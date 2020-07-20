c=========================================================================
c
c     This code solves the QG PV equation in two layers of thicknesses H1 and H2,
c     with the uniform mean background current of speed U_0 in the top layer.
c     Topography can be easily added, as long as it does not violate the QG
c     assumptions (see the code).
c
c     The code requires initial conditions and a set parameters (listed below).
c     The solution is in the form of streamfunction PSI in two layers.
c     At the end of the execution, the program saves the time-mean streamfunction
c     and the restart file "qg2_restart.dta". During teh execution, the program
c     will save "snapshots" (instantaneous fields) at specified intervals. Other
c     diagnostics can be added as well.
c
c     To compile, use, for example,  "ifort qg_channel_students.F -o run.out"
c
c     Author: P. Berloff, Imperial College (London, UK)
c
c=========================================================================
      implicit none

      CHARACTER*2 OTN1
      CHARACTER*3 OTN2
      CHARACTER*4 OTN3
      CHARACTER*30 PSIFILE1
      CHARACTER*30 PSIFILE2


      integer ips,max_time,k_out,k_save,istart,max_spin,ii,jj,ii1,jj1
     & ,ii2,jj2,islip,i,j,l,l_day,m,kk,k_o,k_s,k_day,max_spin1,iic,jjc
     & ,max_time1,l_tot,l_tot_day,leapfrog,i1,ngpswk,idum,k_lagr, k_l
     $     ,ncoarse, max_cycle, ncycle, nrec, n, npoints, iav

      real*8 basinscale,H1,H2,Rd,visc,visc_bot,U_0,h_bot
     & ,weight1,weight2,dt
c
c--- IMPORTANT PARAMETERS (THAT YOU CAN CHANGE):
C              MAX_TIME --- length of integration in days
c              BASINSCALE --- meridional width of the channel (cm)
c              H1,H2 --- layer thickness (cm)
c              Rd  --- Rossby deformation radius (cm)
c              VISC --- lateral Laplacian viscosity (cm**2/sec)
c              VISC_BOT --- bottom friction
c              U_0 --- speed of the background current in the top layer (cm/sec)
c              ISTART --- restart flag (=0 the code starts from analytical IC, 
c                                       =1 the code continues from "qg2_restart.dta" 
c              DT --- time step (in sec). Try making it smaller if the code blows up
c
      parameter(ips=1,max_time=1200,max_spin=0, max_cycle=1, iav=1
     & ,ii=512,jj=256+1)
      parameter(ncoarse=1)
      parameter(basinscale=1800.D5
     & ,H1=1.D5,H2=3.D5,Rd=25.D5
     & ,visc=50.D4
     & ,visc_bot=1.D-7  ! for 4.D-7 spin-down time: 29 days
     & ,U_0=6.D0
     & ,h_bot=.5D0
     & ,istart=0
     & ,islip=1)

      data dt /1080./

      parameter(ii1=ii-1,jj1=jj-1,ii2=ii-2,jj2=jj-2
     & ,weight1=0.05D0,weight2=0.90D0)

c
c--- TIME COUNTERS DETERMINE HOW OFTEN (IN DAYS) VARIOUS OUTPUT IS SAVED
c
      parameter(k_out=10)      ! save data every k_out days
      parameter(k_save=100)    ! save restart every k_save days
c=====================================================================
c
c--- VARIABLES ARE DECLARED BELOW
c

      real*8 psi1(ii,jj,2),psi2(ii,jj,2)
     & ,phi1(ii,jj,2),phi2(ii,jj,2)
     & ,rel1(ii,jj,2),rel2(ii,jj,2)

     & ,zeta1(ii,jj,2),zeta2(ii,jj,2)
     & ,zeta1old(ii,jj,2),zeta2old(ii,jj,2)
     & ,z1(ii,jj,2),z2(ii,jj,2)
     & ,z1old(ii,jj,2),z2old(ii,jj,2)

     & ,pv1(0:ii+1,jj),pv2(0:ii+1,jj),pv1_x(ii,jj)
     & ,q1(0:ii+1,jj),q2(0:ii+1,jj)
     & ,u1(0:ii+1,jj),u2(0:ii+1,jj)
     & ,v1(0:ii+1,jj),v2(0:ii+1,jj)
     & ,diss1(ii,jj),diss2(ii,jj)

     & ,b1(jj2,ii),b2(jj2,ii)
     & ,a_ell(jj2),d_ell(jj2),u_ell(jj2)
     & ,phi1_corr1(ii,jj),phi2_corr1(ii,jj)
     & ,phi1_corr2(ii,jj),phi2_corr2(ii,jj)
     & ,z1_ell(ii,jj),z2_ell(ii,jj)

     & ,beta_y(jj),beta1_y(jj),beta2_y(jj)
     & ,bot(0:ii+1,jj),bot_x(ii-1,jj-1),bot_y(ii-1,jj-1)
     & ,S1,S2,SS,f0,hh,TT
     & ,theta11,theta12,theta21,theta22
     & ,omega11,omega12,omega21,omega22
     & ,rho1_11,rho1_12,rho1_21,rho1_22
     & ,rho2_11,rho2_12,rho2_21,rho2_22
     & ,kappa11,kappa12,kappa21,kappa22
     & ,D1_s,D2_s,D1_n,D2_n,A1_s,A2_s,A1_n,A2_n
     & ,C1_s,C2_s,C1_n,C2_n
     & ,alpha1,alpha2,beta1,beta2
     & ,alpha1_s,alpha2_s,alpha1_n,alpha2_n
     & ,beta1_s,beta2_s,beta1_n,beta2_n
     & ,pi,dt_day,uscale,scale,tscale,beta,cff,cff1,cff2
     & ,beta_nondim,visc_nondim,visc_bot_nondim
     & ,time_day,time_day1,time0_day
     & ,ekin,ekin1,ekin2,epot
     & ,z1_av,z2_av, dt2,dt6
     $ ,psi1_av(0:ii+2,jj), psi2_av(0:ii+2,jj) 
      real aux(ii,jj),aux1(ii,jj),aux2(ii,jj),aux3(ii,jj)
     &               ,aux4(ii,jj),aux5(ii,jj),aux6(ii,jj)
     & ,unused,ran1


c
c     Diagnistics (r4, for Matlab)
c
      real*4 psi1s(ii,jj), psi2s(ii,jj) 
      external ran1
C============================================================

      data pi/3.14159265358979323846D0/
      data beta/2.D-13/

      common /ONE/ psi1,psi2
      common /TWO/ z1,z2,phi1,phi2
      common /THREE/ rel1,rel2
      common /FOUR/ theta11,theta12,theta21,theta22
      common /FIVE/ omega11,omega12,omega21,omega22
      common /SIX/ rho1_11,rho1_12,rho1_21,rho1_22
     &            ,rho2_11,rho2_12,rho2_21,rho2_22
      common /SEVEN/ kappa11,kappa12,kappa21,kappa22
c
c--- nondimensionalization
c
      uscale=1.
      scale=basinscale/dfloat(jj-1)
      tscale=scale/uscale

      dt=dt/tscale
      dt2=dt/2.0
      dt6=dt/6.0
      dt_day=86400./tscale

      beta_nondim=beta*scale*scale/uscale
      visc_nondim=visc/(scale*uscale)
      visc_bot_nondim=visc_bot*scale/uscale
c
c--- stratification
c
      SS=(scale/Rd)**2

      S1=SS/(1.+H1/H2)
      S2=SS-S1
      write(*,*)'S1,S2=',S1,S2

      k_o=0
      k_s=0
      k_l=0

      TT=2.*dt

      l_tot=int(86400./(dt*tscale)+0.001)
      l_tot_day=int(86400./(dt*tscale)+0.001)
      l=0
      l_day=0

      max_time1=max_time*l_tot_day/2
      max_spin1=max_spin*l_tot_day/2

      time_day=0.
c
c--- TOPOGRAPHY CAN BE PRESCRIBED HERE
c
      f0=.83D-4
      hh=h_bot*f0*tscale

      do j=1,jj
         do i=1,ii
            bot(i,j)=0.
         enddo
      enddo

      do j=1,jj
         do i=ii/4+1,2*ii/4
c           bot(i,j)=hh*dsin(2.D0*pi*float(i-ii/4-1)/dfloat(ii/2))
         enddo
      enddo

c      do j=1,jj
c         do i=1,ii
c            if(i.eq.1)then
c              bot_x(i,j)=.5*(bot(i+1,j)-bot(ii,j))
c            elseif(i.eq.ii)then
c              bot_x(i,j)=.5*(bot(1,j)-bot(i-1,j))
c            else
c              bot_x(i,j)=.5*(bot(i+1,j)-bot(i-1,j))
c            endif
c         enddo
c      enddo
c
c      do i=1,ii
c         do j=1,jj
c            if(j.eq.1)then
c              bot_y(i,j)=bot(i,j+1)-bot(i,j)
c            elseif(j.eq.jj)then
c              bot_y(i,j)=bot(i,j)-bot(i,j-1)
c            else
c              bot_y(i,j)=.5*(bot(i,j+1)-bot(i,j-1))
c            endif
c         enddo
c      enddo
c
c--- beta*y
c
      do j=1,jj
         beta_y(j)=beta_nondim*(dfloat(j-1)-dfloat(jj-1)/2.D0)
         beta1_y(j)=(beta_nondim+S1*U_0)*(dfloat(j-1)-dfloat(jj-1)/2.D0)
         beta2_y(j)=(beta_nondim-S2*U_0)*(dfloat(j-1)-dfloat(jj-1)/2.D0)
      enddo
c
c--- initialization for elliptic solver
c
      theta11=H1/(H1+H2)
      theta12=H2/(H1+H2)
      theta21= 1.
      theta22=-1.

      omega11=1.
      omega12= H2/(H1+H2)
      omega21=1.
      omega22=-H1/(H1+H2)
c
c--- initialization of the homogeneous elliptic problem
c
      do j=1,jj2
         a_ell(j)=1.
         b1(j,1)=-2.D0
         b2(j,1)=-2.D0-SS
      enddo

      do i=1,ii/2
         do j=1,jj2
           i1=2*i-1
           b1(j,i1)=-2.*(2.-dcos(2.D0*pi*(i-1)/ii))
           b2(j,i1)=-SS+b1(j,i1)
           b1(j,i1+1)=b1(j,i1)
           b2(j,i1+1)=b2(j,i1)
         enddo
         do j=1,jj2
            b1(j,2)=-2.*(2.-dcos(2.D0*pi*(ii/2+1-1)/ii))
            b2(j,2)=-SS+b1(j,2)
         enddo 
      enddo
c
c--- inhomogeneous elliptic problems with (1,0) and (0,1) BCs
c
      do j=2,jj2
         d_ell(j)=0.D0
      enddo
      d_ell(1)=-1.
      call tridag(a_ell,b1(1,1),a_ell,d_ell,u_ell,jj2)
      do i=1,ii
         do j=2,jj-1
            phi1_corr1(i,j)=u_ell(j-1)
         enddo
         phi1_corr1(i,1)=1.D0
         phi1_corr1(i,jj)=0.D0
      enddo

      do j=1,jj2-1
         d_ell(j)=0.D0
      enddo
      d_ell(jj2)=-1.
      call tridag(a_ell,b1(1,1),a_ell,d_ell,u_ell,jj2)
      do i=1,ii
         do j=2,jj-1
            phi1_corr2(i,j)=u_ell(j-1)
         enddo
         phi1_corr2(i,1)=0.D0
         phi1_corr2(i,jj)=1.D0
      enddo
c
c---
c
      do j=2,jj2
         d_ell(j)=0.D0
      enddo
      d_ell(1)=-1.
      call tridag(a_ell,b2(1,1),a_ell,d_ell,u_ell,jj2)
      do i=1,ii
         do j=2,jj-1
            phi2_corr1(i,j)=u_ell(j-1)
         enddo
         phi2_corr1(i,1)=1.D0
         phi2_corr1(i,jj)=0.D0
      enddo

      do j=1,jj2-1
         d_ell(j)=0.D0
      enddo
      d_ell(jj2)=-1.
      call tridag(a_ell,b2(1,1),a_ell,d_ell,u_ell,jj2)
      do i=1,ii
         do j=2,jj-1
            phi2_corr2(i,j)=u_ell(j-1)
         enddo
         phi2_corr2(i,1)=0.D0
         phi2_corr2(i,jj)=1.D0
      enddo
c
c--- viscous boundary fluxes of the corrections
c
      C1_s=0.
      C2_s=0.
      C1_n=0.
      C2_n=0.
      do i=1,ii
         C1_s=C1_s-3.*phi1_corr1(i,1)+4.*phi1_corr1(i,2)
     &                                -phi1_corr1(i,3)
         C2_s=C2_s-3.*phi2_corr1(i,1)+4.*phi2_corr1(i,2)
     &                                -phi2_corr1(i,3)
         C1_n=C1_n-3.*phi1_corr1(i,jj)+4.*phi1_corr1(i,jj-1)
     &                                -phi1_corr1(i,jj-2)
         C2_n=C2_n-3.*phi2_corr1(i,jj)+4.*phi2_corr1(i,jj-1)
     &                                -phi2_corr1(i,jj-2)
      enddo
      C1_s= visc_nondim*C1_s
      C2_s= visc_nondim*C2_s
      C1_n= visc_nondim*C1_n
      C2_n= visc_nondim*C2_n
      write(*,*)'C1_s=',C1_s,'; C1_n=',C1_n
      write(*,*)'C2_s=',C2_s,'; C2_n=',C2_n

      alpha1_s=-1./C1_s  ! analytically, all C_s,n have to be ZEROS as the
      alpha1_n=-1./C1_n  ! third derivatives of a linear function (and no
      beta1_n=-1./C1_s   ! alpha and beta corrections were needed), but
      beta1_s=-1./C1_n   ! numerically they are not zeros!

      alpha2_s=-1./C2_s
      alpha2_n=0.
      beta2_s=0.
      beta2_n=-1./C2_s

      write(*,*)'alpha1_s,n=',alpha1_s,alpha1_n
      write(*,*)'beta1_s,n=',beta1_s,beta1_n
      write(*,*)'alpha2_s,n=',alpha2_s,alpha2_n
      write(*,*)'beta2_s,n=',beta2_s,beta2_n
c
c--- Jacobian-coefficients of the equations written with the vertical modes
c
      rho1_11=theta11*omega11*omega11+theta12*omega21*omega21
      rho1_12=theta11*omega11*omega12+theta12*omega21*omega22

      rho1_21=theta11*omega12*omega11+theta12*omega22*omega21
      rho1_22=theta11*omega12*omega12+theta12*omega22*omega22

      write(*,*)'rho1_11 _12=',rho1_11,rho1_12
      write(*,*)'rho1_21 _22=',rho1_21,rho1_22
      write(*,*)'                                          '

      rho2_11=theta21*omega11*omega11+theta22*omega21*omega21
      rho2_12=theta21*omega11*omega12+theta22*omega21*omega22

      rho2_21=theta21*omega12*omega11+theta22*omega22*omega21
      rho2_22=theta21*omega12*omega12+theta22*omega22*omega22

      write(*,*)'rho2_11 _12=',rho2_11,rho2_12
      write(*,*)'rho2_21 _22=',rho2_21,rho2_22
      write(*,*)'                                          '
c
c--- bottom friction coefficients for vertical modes
c
c      kappa11=visc_bot_nondim*theta12*omega21
c      kappa12=visc_bot_nondim*theta12*omega22
c      kappa21=visc_bot_nondim*theta22*omega21
c      kappa22=visc_bot_nondim*theta22*omega22
c      write(*,*)'kappa11=',kappa11,'; kappa12=',kappa12
c      write(*,*)'kappa21=',kappa21,'; kappa22=',kappa22
c
c--- initial conditions
c

       if(istart.eq.0)then
        do j=1,jj
           do i=1,ii
c              psi1(i,j,1)=0.
c              psi1(i,j,2)=0.
c              psi2(i,j,1)=0.
c              psi2(i,j,2)=0.

              psi1(i,j,1)=10.*dsin(2.*pi*7.7*float(i-1)/float(jj))
              psi2(i,j,1)=.25*psi1(i,j,1)  ! U=+6,0

c              psi1(i,j,1)=3.*dsin(2.*pi*8.9*float(i-1)/float(jj))
c              psi2(i,j,1)=-psi1(i,j,1)/.3  ! U=-3,0

c              psi1(i,j,1)=4.*ran1(idum)

              psi1(i,j,2)=psi1(i,j,1)
              psi2(i,j,2)=psi2(i,j,1)

c              write(*,*)'Test',psi1(i,j,1),psi2(i,j,1)

           enddo
        enddo

        call subtract_mean_periodic(ii,jj,psi1)
        call subtract_mean_periodic(ii,jj,psi2)

      else
        open(55,file='qg2_restart.dta',form='unformatted')
        read(55) psi1,psi2
        close(55)
      endif

      call rel_from_psi(islip,ii,jj,psi1,psi2,rel1,rel2)
      call zeta_from_psi_and_rel(ii,jj,S1,S2,psi1,psi2,rel1,rel2
     &                                                ,zeta1,zeta2)

      do j=1,jj
         do i=1,ii
            zeta1old(i,j,1)=zeta1(i,j,1)
            zeta2old(i,j,1)=zeta2(i,j,1)
         enddo
      enddo
c
c--- zeros
c
      do j=1,jj
         do i=1,ii
            psi1_av(i,j)=0.
            psi2_av(i,j)=0.
         enddo
      enddo

      do j=1,jj
         do i=0,ii+1
            u1(i,j)=0.
            u2(i,j)=0.
            v1(i,j)=0.
            v2(i,j)=0.
         enddo
      enddo


c
c-- Opening eddy energy files
c
      open(unit=11,file='outputs/eddy_energy.dat',status='unknown')




c
cccccccccccccccccccccccccccccccccccccccccccccc
      do 2000 ncycle=1,max_cycle
         write(*,*)'############## CYCLE:',ncycle
         k_o=0
         k_l=0
         l=0
         l_day=0
         time_day=0.
         nrec=0

c
c################################# MAIN CYCLE #################################
c
      do 1000 kk=1,max_spin1+max_time1
      do 999 leapfrog=1,2
c         write(*,*)'kk=',kk,'; leapfrog=',leapfrog
         time_day=time_day+dt*tscale/86400.
         l=l+1
         l_day=l_day+1

         k_day=int(time_day+0.001)

         if(leapfrog.eq.1)then
           m=2
         elseif(leapfrog.eq.2)then
           m=1
         endif

         if(l_day.eq.l_tot_day)then
           l_day=0
           k_o=k_o+1
           k_s=k_s+1
           k_l=k_l+1
         endif

         if(l.eq.l_tot)then
           l=0
         endif
c
c--- zeta in the interior
c
         call zeta_interior(leapfrog,m,TT,ii,jj,S1,S2,U_0
     &                     ,beta1_y,beta2_y
     &                     ,bot
     &                     ,visc_nondim,visc_bot_nondim
     &                     ,psi1,psi2,rel1,rel2
     &                     ,pv1,pv2,pv1_x,q1,q2,u1,u2,v1,v2
     &                     ,diss1,diss2
     &                     ,zeta1,zeta2
     &                     ,D1_s,D2_s,D1_n,D2_n
     &                     ,A1_s,A2_s,A1_n,A2_n)
c
c--- solve helmholtz problems
c
         do j=1,jj
            do i=1,ii
               z1_ell(i,j)= theta11*zeta1(i,j,leapfrog)
     &                     +theta12*zeta2(i,j,leapfrog)
               z2_ell(i,j)= theta21*zeta1(i,j,leapfrog)
     &                     +theta22*zeta2(i,j,leapfrog)
            enddo
         enddo
         do i=1,ii  ! needed by the elliptic solver
            phi1(i,1,leapfrog)=0.
            phi2(i,1,leapfrog)=0.
            phi1(i,jj,leapfrog)=0.
            phi2(i,jj,leapfrog)=0.
         enddo

         call chsolv(ii,jj,z1_ell,phi1(1,1,leapfrog)
     &                                     ,a_ell,b1,d_ell,u_ell)
         call chsolv(ii,jj,z2_ell,phi2(1,1,leapfrog)
     &                                     ,a_ell,b2,d_ell,u_ell)
c
c--- convert to layers
c
         do j=1,jj
            do i=1,ii
               psi1(i,j,leapfrog)= omega11*phi1(i,j,leapfrog)
     &                            +omega12*phi2(i,j,leapfrog)
               psi2(i,j,leapfrog)= omega21*phi1(i,j,leapfrog)
     &                            +omega22*phi2(i,j,leapfrog)
            enddo
         enddo
c
c--- enforce momentum constraint
c
         call phi_bry_fluxes(leapfrog,m,ii,jj,S1,S2
     &                            ,beta1_y,beta2_y,visc_nondim
     &                            ,psi1,psi2
     &                            ,pv1,pv2
     &                            ,q1,q2
     &                            ,v1,v2
     &                            ,A1_s,A2_s,A1_n,A2_n
     &                            ,D1_s,D2_s,D1_n,D2_n
     &                            ,C1_s,C2_s,C1_n,C2_n)

         alpha1=alpha1_s*C1_s+alpha1_n*C1_n
         alpha2=alpha2_s*C2_s+alpha2_n*C2_n
         beta1=beta1_s*C1_s+beta1_n*C1_n
         beta2=beta2_s*C2_s+beta2_n*C2_n

         do j=1,jj
            do i=1,ii
               phi2(i,j,leapfrog)=phi2(i,j,leapfrog)
     &               +alpha2*phi2_corr1(i,j)+beta2*phi2_corr2(i,j)
            enddo
         enddo
c
c--- enforce mass constraint
c
         call subtract_mean_periodic(ii,jj,phi2(1,1,leapfrog))
c
c--- AGAIN convert to layers
c
         do j=1,jj
            do i=1,ii
               psi1(i,j,leapfrog)= omega11*phi1(i,j,leapfrog)
     &                            +omega12*phi2(i,j,leapfrog)
               psi2(i,j,leapfrog)= omega21*phi1(i,j,leapfrog)
     &                            +omega22*phi2(i,j,leapfrog)
            enddo
         enddo
c
c--- smooth z in time
c
         if(leapfrog.eq.1)then
           do j=1,jj
              do i=1,ii
                 zeta1(i,j,2)=weight1*(zeta1old(i,j,1)+zeta1(i,j,1))
     &                    +weight2*zeta1(i,j,2)
                 zeta2(i,j,2)=weight1*(zeta2old(i,j,1)+zeta2(i,j,1))
     &                    +weight2*zeta2(i,j,2)
                 zeta1old(i,j,2)=zeta1(i,j,2)
                 zeta2old(i,j,2)=zeta2(i,j,2)
              enddo
           enddo
         endif
         if(leapfrog.eq.2)then
           do j=1,jj
              do i=1,ii
                 zeta1(i,j,1)=weight1*(zeta1old(i,j,2)+zeta1(i,j,2))
     &                    +weight2*zeta1(i,j,1)
                 zeta2(i,j,1)=weight1*(zeta2old(i,j,2)+zeta2(i,j,2))
     &                    +weight2*zeta2(i,j,1)
                 zeta1old(i,j,1)=zeta1(i,j,1)
                 zeta2old(i,j,1)=zeta2(i,j,1)
              enddo
           enddo
         endif
c
c--- relative vorticity
c
         do j=1,jj
            do i=1,ii
               rel1(i,j,leapfrog)=zeta1(i,j,leapfrog)
     &                  +S1*(psi1(i,j,leapfrog)-psi2(i,j,leapfrog))
               rel2(i,j,leapfrog)=zeta2(i,j,leapfrog)
     &                  +S2*(psi2(i,j,leapfrog)-psi1(i,j,leapfrog))
            enddo
         enddo
         call boundary_condition(islip,leapfrog,ii,jj
     &                                ,rel1,rel2,psi1,psi2)
c
c--- time-mean streamfunction
c
         if(kk.gt.max_spin1)then
           do j=1,jj
              do i=1,ii
                 psi1_av(i,j)=psi1_av(i,j)+psi1(i,j,leapfrog)
                 psi2_av(i,j)=psi2_av(i,j)+psi2(i,j,leapfrog)
              enddo
           enddo
         endif
c
c--- diagnostics
c
         if(k_o.eq.k_out)then
           k_o=0
           if(kk.le.max_spin1)then
             write(*,*)'Spin-up: Time (days) =',k_day
           else
             write(*,*)'Time (days) =',k_day
           endif

           call energy(H1,H2,S1,S2,ii,jj,psi1,psi2,ekin1,ekin2,epot)

           write(*,'(A5,F12.6,A8,F12.6,A8,F12.6)')
     & 'epot=',epot,'; ekin1=',ekin1,'; ekin2=',ekin2


            write(11,*) epot,ekin1,ekin2


c
c-- Saving Snapshots
c

            if(k_day.lt.100)then
                WRITE(OTN1,666),k_day
            else if (k_day.ge.1000)then
                WRITE(OTN3,888),k_day
            else
                WRITE(OTN2,777),k_day
            endif


666         FORMAT(I2)

777         FORMAT(I3)

888         FORMAT(I4)


            if(k_day.lt.100)then
                PSIFILE1='outputs/psi1/psi1_0'//OTN1//'.bin'
                PSIFILE2='outputs/psi2/psi2_0'//OTN1//'.bin'
            else if (k_day.ge.1000)then
                PSIFILE1='outputs/psi1/psi1_'//OTN3//'.bin'
                PSIFILE2='outputs/psi2/psi2_'//OTN3//'.bin'
            else
                PSIFILE1='outputs/psi1/psi1_'//OTN2//'.bin'
                PSIFILE2='outputs/psi2/psi2_'//OTN2//'.bin'
            endif

            open(31,file=PSIFILE1,form='unformatted')
            open(32,file=PSIFILE2,form='unformatted')


            psi1s=psi1(:,:,2)
            psi2s=psi2(:,:,2)


               write(31) psi1s
               write(32) psi2s


            close(31)
            close(32)

         endif
c
c-- save the restart file
c
         if(k_s.eq.k_save)then
           k_s=0

           open(55,file='qg2_restart.dta',form='unformatted')
           write(55) psi1,psi2
           close(55)

           write(*,*)'Output plotted...'
        endif


 999  continue
 1000 continue
 2000 continue
c
c--- time-mean streamfunction
c
      do j=1,jj
         do i=0,ii+2
            psi1_av(i,j)=psi1_av(i,j)/float(2*max_time1*max_cycle)
            psi2_av(i,j)=psi2_av(i,j)/float(2*max_time1*max_cycle)
         enddo
      enddo

      open(41,file='outputs/qg2_aver.dat',form='unformatted')
      write(41) psi1_av,psi2_av
      close(41)


       close(11)


      stop
      end

C============================================================================
C
      subroutine boundary_condition(islip,leapfrog,ii,jj
     &                                ,rel1,rel2,psi1,psi2)
      implicit none
      integer ii,jj,i,j,islip,leapfrog
      real*8 rel1(ii,jj,2),rel2(ii,jj,2)
     &      ,psi1(ii,jj,2),psi2(ii,jj,2)
       
      if(islip.eq.1)then         ! NO-SLIP
        do i=1,ii
c           rel1(i, 1,leapfrog)=-3.5*psi1(i,1   ,leapfrog)
c     &                         +4.0*psi1(i,2   ,leapfrog)
c     &                         -0.5*psi1(i,3   ,leapfrog)
c           rel1(i,jj,leapfrog)=-3.5*psi1(i,jj  ,leapfrog)
c     &                         +4.0*psi1(i,jj-1,leapfrog)
c     &                         -0.5*psi1(i,jj-2,leapfrog)
c
c           rel2(i, 1,leapfrog)=-3.5*psi2(i,1   ,leapfrog)
c     &                         +4.0*psi2(i,2   ,leapfrog)
c     &                         -0.5*psi2(i,3   ,leapfrog)
c           rel2(i,jj,leapfrog)=-3.5*psi2(i,jj  ,leapfrog)
c     &                         +4.0*psi2(i,jj-1,leapfrog)
c     &                         -0.5*psi2(i,jj-2,leapfrog)

           rel1(i, 1,leapfrog)=2.*(-psi1(i,1   ,leapfrog)
     &                             +psi1(i,2   ,leapfrog))
           rel1(i,jj,leapfrog)=2.*(-psi1(i,jj  ,leapfrog)
     &                             +psi1(i,jj-1,leapfrog))

           rel2(i, 1,leapfrog)=2.*(-psi2(i,1   ,leapfrog)
     &                             +psi2(i,2   ,leapfrog))
           rel2(i,jj,leapfrog)=2.*(-psi2(i,jj  ,leapfrog)
     &                             +psi2(i,jj-1,leapfrog))
        enddo
      else                       ! FREE-SLIP
        do i=1,ii
           rel1(i, 1,leapfrog)=0.
           rel1(i,jj,leapfrog)=0.

           rel2(i, 1,leapfrog)=0.
           rel2(i,jj,leapfrog)=0.
        enddo
      endif

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine phi_bry_fluxes(leapfrog,m,ii,jj,S1,S2
     &                            ,beta1_y,beta2_y,visc
     &                            ,psi1,psi2
     &                            ,pv1,pv2
     &                            ,q1,q2
     &                            ,v1,v2
     &                            ,A1_s,A2_s,A1_n,A2_n
     &                            ,D1_s,D2_s,D1_n,D2_n
     &                            ,C1_s,C2_s,C1_n,C2_n)
      implicit none
      integer leapfrog,m,ii,jj,i,j
      real*8 psi1(ii,jj,2),psi2(ii,jj,2)
     & ,beta1_y(jj),beta2_y(jj),visc

     & ,pv1(0:ii+1,jj),pv2(0:ii+1,jj)
     & ,q1(0:ii+1,jj),q2(0:ii+1,jj)
     & ,v1(0:ii+1,jj),v2(0:ii+1,jj)

     & ,S1,S2
     & ,theta11,theta12,theta21,theta22
     & ,rho1_11,rho1_12,rho1_21,rho1_22
     & ,rho2_11,rho2_12,rho2_21,rho2_22
     & ,D1_s,D2_s,D1_n,D2_n
     & ,A1_s,A2_s,A1_n,A2_n
     & ,C1_s,C2_s,C1_n,C2_n

      common /FOUR/ theta11,theta12,theta21,theta22
      common /SIX/ rho1_11,rho1_12,rho1_21,rho1_22
     &            ,rho2_11,rho2_12,rho2_21,rho2_22
c
c--- auxiliary arrays
c
      do i=2,ii-1
         pv1(i,2)=-4.*psi1(i,2,leapfrog)
     &              +psi1(i-1,2,leapfrog)+psi1(i+1,2,leapfrog)
     &              +psi1(i,1,leapfrog)+psi1(i,3,leapfrog)
     &            -S1*(psi1(i,2,leapfrog)-psi2(i,2,leapfrog))
         pv2(i,2)=-4.*psi2(i,2,leapfrog)
     &              +psi2(i-1,2,leapfrog)+psi2(i+1,2,leapfrog)
     &              +psi2(i,1,leapfrog)+psi2(i,3,leapfrog)
     &            -S2*(psi2(i,2,leapfrog)-psi1(i,2,leapfrog))

         pv1(i,jj-1)=-4.*psi1(i,jj-1,leapfrog)
     &              +psi1(i-1,jj-1,leapfrog)+psi1(i+1,jj-1,leapfrog)
     &              +psi1(i,jj-2,leapfrog)+psi1(i,jj,leapfrog)
     &            -S1*(psi1(i,jj-1,leapfrog)-psi2(i,jj-1,leapfrog))
         pv2(i,jj-1)=-4.*psi2(i,jj-1,leapfrog)
     &              +psi2(i-1,jj-1,leapfrog)+psi2(i+1,jj-1,leapfrog)
     &              +psi2(i,jj-2,leapfrog)+psi2(i,jj,leapfrog)
     &            -S2*(psi2(i,jj-1,leapfrog)-psi1(i,jj-1,leapfrog))
      enddo

      j=2
      pv1(1,j)=-4.*psi1(1,j,leapfrog)
     &              +psi1(ii,j,leapfrog)+psi1(2,j,leapfrog)
     &              +psi1(1,j-1,leapfrog)+psi1(1,j+1,leapfrog)
     &            -S1*(psi1(1,j,leapfrog)-psi2(1,j,leapfrog))
      pv2(1,j)=-4.*psi2(1,j,leapfrog)
     &              +psi2(ii,j,leapfrog)+psi2(2,j,leapfrog)
     &              +psi2(1,j-1,leapfrog)+psi2(1,j+1,leapfrog)
     &            -S2*(psi2(1,j,leapfrog)-psi1(1,j,leapfrog))
      pv1(ii,j)=-4.*psi1(ii,j,leapfrog)
     &              +psi1(ii-1,j,leapfrog)+psi1(1,j,leapfrog)
     &              +psi1(ii,j-1,leapfrog)+psi1(ii,j+1,leapfrog)
     &            -S1*(psi1(ii,j,leapfrog)-psi2(ii,j,leapfrog))
      pv2(ii,j)=-4.*psi2(ii,j,leapfrog)
     &              +psi2(ii-1,j,leapfrog)+psi2(1,j,leapfrog)
     &              +psi2(ii,j-1,leapfrog)+psi2(ii,j+1,leapfrog)
     &            -S2*(psi2(ii,j,leapfrog)-psi1(ii,j,leapfrog))

      j=jj-1
      pv1(1,j)=-4.*psi1(1,j,leapfrog)
     &              +psi1(ii,j,leapfrog)+psi1(2,j,leapfrog)
     &              +psi1(1,j-1,leapfrog)+psi1(1,j+1,leapfrog)
     &            -S1*(psi1(1,j,leapfrog)-psi2(1,j,leapfrog))
      pv2(1,j)=-4.*psi2(1,j,leapfrog)
     &              +psi2(ii,j,leapfrog)+psi2(2,j,leapfrog)
     &              +psi2(1,j-1,leapfrog)+psi2(1,j+1,leapfrog)
     &            -S2*(psi2(1,j,leapfrog)-psi1(1,j,leapfrog))
      pv1(ii,j)=-4.*psi1(ii,j,leapfrog)
     &              +psi1(ii-1,j,leapfrog)+psi1(1,j,leapfrog)
     &              +psi1(ii,j-1,leapfrog)+psi1(ii,j+1,leapfrog)
     &            -S1*(psi1(ii,j,leapfrog)-psi2(ii,j,leapfrog))
      pv2(ii,j)=-4.*psi2(ii,j,leapfrog)
     &              +psi2(ii-1,j,leapfrog)+psi2(1,j,leapfrog)
     &              +psi2(ii,j-1,leapfrog)+psi2(ii,j+1,leapfrog)
     &            -S2*(psi2(ii,j,leapfrog)-psi1(ii,j,leapfrog))

      do i=0,ii+1
         q1(i,2)=pv1(i,2)+beta1_y(2)
         q2(i,2)=pv2(i,2)+beta2_y(2)
         q1(i,jj-1)=pv1(i,jj-1)+beta1_y(jj-1)
         q2(i,jj-1)=pv2(i,jj-1)+beta2_y(jj-1)
      enddo

      do i=2,ii-1
         v1(i,2)=psi1(i+1,2,leapfrog)-psi1(i-1,2,leapfrog)
         v2(i,2)=psi2(i+1,2,leapfrog)-psi2(i-1,2,leapfrog)
         v1(i,jj-1)=psi1(i+1,jj-1,leapfrog)-psi1(i-1,jj-1,leapfrog)
         v2(i,jj-1)=psi2(i+1,jj-1,leapfrog)-psi2(i-1,jj-1,leapfrog)
      enddo
      v1(1,2)=psi1(2,2,leapfrog)-psi1(ii,2,leapfrog)
      v2(1,2)=psi2(2,2,leapfrog)-psi2(ii,2,leapfrog)
      v1(ii,2)=psi1(1,2,leapfrog)-psi1(ii-1,2,leapfrog)
      v2(ii,2)=psi2(1,2,leapfrog)-psi2(ii-1,2,leapfrog)
      v1(1,jj-1)=psi1(2,jj-1,leapfrog)-psi1(ii,jj-1,leapfrog)
      v2(1,jj-1)=psi2(2,jj-1,leapfrog)-psi2(ii,jj-1,leapfrog)
      v1(ii,jj-1)=psi1(1,jj-1,leapfrog)-psi1(ii-1,jj-1,leapfrog)
      v2(ii,jj-1)=psi2(1,jj-1,leapfrog)-psi2(ii-1,jj-1,leapfrog)
c
c--- viscous boundary fluxes
c
      D1_s=0.
      D2_s=0.
      D1_n=0.
      D2_n=0.
      do i=1,ii
         D1_s=D1_s+visc*(
     &             -3.*(theta11*psi1(i,1,m)+theta12*psi2(i,1,m))
     &             +4.*(theta11*psi1(i,2,m)+theta12*psi2(i,2,m))
     &                -(theta11*psi1(i,3,m)+theta12*psi2(i,3,m))
     &                  )
         D2_s=D2_s+visc*(
     &             -3.*(theta21*psi1(i,1,m)+theta22*psi2(i,1,m))
     &             +4.*(theta21*psi1(i,2,m)+theta22*psi2(i,2,m))
     &                -(theta21*psi1(i,3,m)+theta22*psi2(i,3,m))
     &                  )
         D1_n=D1_n+visc*(
     &             -3.*(theta11*psi1(i,jj,m)+theta12*psi2(i,jj,m))
     &             +4.*(theta11*psi1(i,jj-1,m)+theta12*psi2(i,jj-1,m))
     &                -(theta11*psi1(i,jj-2,m)+theta12*psi2(i,jj-2,m))
     &                  )
         D2_n=D2_n+visc*(
     &             -3.*(theta21*psi1(i,jj,m)+theta22*psi2(i,jj,m))
     &             +4.*(theta21*psi1(i,jj-1,m)+theta22*psi2(i,jj-1,m))
     &                -(theta21*psi1(i,jj-2,m)+theta22*psi2(i,jj-2,m))
     &                  )
      enddo
c
c--- advective boundary fluxes
c
      A1_s=0.
      A2_s=0.
      A1_n=0.
      A2_n=0.
      do i=1,ii
         A1_s=A1_s+rho1_11*(theta11*v1(i,2)+theta12*v2(i,2))
     &                    *(theta11*q1(i,2)+theta12*q2(i,2))
     &            +rho1_22*(theta21*v1(i,2)+theta22*v2(i,2))
     &                    *(theta21*pv1(i,2)+theta22*pv2(i,2))
         A1_n=A1_n-rho1_11*(theta11*v1(i,jj-1)+theta12*v2(i,jj-1))
     &                    *(theta11*q1(i,jj-1)+theta12*q2(i,jj-1))
     &            -rho1_22*(theta21*v1(i,jj-1)+theta22*v2(i,jj-1))
     &                    *(theta21*pv1(i,jj-1)+theta22*pv2(i,jj-1))

         A2_s=A2_s+rho2_12*(theta11*v1(i,2)+theta12*v2(i,2))
     &                    *(theta21*pv1(i,2)+theta22*pv2(i,2))
     &            +rho2_21*(theta21*v1(i,2)+theta22*v2(i,2))
     &                    *(theta11*pv1(i,2)+theta12*pv2(i,2))
     &            +rho2_22*(theta21*v1(i,2)+theta22*v2(i,2))
     &                    *(theta21*q1(i,2)+theta22*q2(i,2))

         A2_n=A2_n-rho2_12*(theta11*v1(i,jj-1)+theta12*v2(i,jj-1))
     &                    *(theta21*pv1(i,jj-1)+theta22*pv2(i,jj-1))
     &            -rho2_21*(theta21*v1(i,jj-1)+theta22*v2(i,jj-1))
     &                    *(theta11*pv1(i,jj-1)+theta12*pv2(i,jj-1))
     &            -rho2_22*(theta21*v1(i,jj-1)+theta22*v2(i,jj-1))
     &                    *(theta21*q1(i,jj-1)+theta22*q2(i,jj-1))
      enddo
      A1_s=.25*A1_s
      A1_n=.25*A1_n
      A2_s=.25*A2_s
      A2_n=.25*A2_n

      C1_s=D1_s+A1_s
      C2_s=D2_s+A2_s
      C1_n=D1_n+A1_n
      C2_n=D2_n+A2_n

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine chsolv(jj,ii,z,phi,a,b,d,u)
      implicit none
      integer ii,jj,i,j
      real*8 z(jj,ii),phi(jj,ii),a(ii-2),b(ii-2,jj),d(ii-2),u(ii-2)

      do i=2,ii-1
         call realft(z(1,i),jj,1)
      enddo

      do j=1,jj
         do i=2,ii-1
            d(i-1)=z(j,i)
         enddo
         call tridag(a,b(1,j),a,d,u,ii-2)
         do i=1,ii-2
            phi(j,i+1)=u(i)
         enddo
      enddo

      do i=2,ii-1
         call realft(phi(1,i),jj,-1)
         do j=1,jj
            phi(j,i)=phi(j,i)*2.D0/float(jj)
         enddo
      enddo

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine realft(data,n,isign)
        implicit real*8 (a-h,o-z)
      INTEGER isign,n
      real*8 data(n)
      INTEGER i,i1,i2,i3,i4,n2p3

      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.1) then
        c2=-0.5
        call four1(data,n/2,+1)
      else
        c2=0.5
        theta=-theta
      endif
      wpr=-2.0d0*dsin(0.5d0*theta)**2
      wpi=dsin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=wr
        wis=wi
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n/2,-1)
      endif

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine four1(data,nn,isign)
        implicit real*8 (a-h,o-z)
      INTEGER isign,nn
      real*8 data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=wr*data(j)-wi*data(j+1)
            tempi=wr*data(j+1)+wi*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine tridag(a,b,c,r,u,n)
        implicit real*8 (a-h,o-z)
      INTEGER NMAX
      PARAMETER (NMAX=1000)
      real*8 a(n),b(n),c(n),r(n),u(n)
c
      INTEGER j
      real*8 gam(NMAX)
      if (n.gt.NMAX-2) PAUSE 'check dimensions'
      if(b(1).eq.0.)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine cubic_time(a,psi1,psi2,psi3,psi4,jj,ii,p)
      implicit none
      integer jj,ii,j,i
      real*8 a,c1,c2,a2,a3,cff1,cff2,cff3,cff4
      real psi1(jj,ii),psi2(jj,ii),psi3(jj,ii),psi4(jj,ii)
     & ,p(jj,ii)
      parameter(c1=1./6.,c2=2./6.)

      a2=a**2
      a3=a**3
      cff1=  -c2*a+.5*a2-c1*a3
      cff2=1.-.5*a   -a2+.5*a3
      cff3=      a+.5*a2-.5*a3
      cff4=  -c1*a      +c1*a3
      do i=1,ii
         do j=1,jj
            p(j,i)=cff1*psi1(j,i)+cff2*psi2(j,i)
     &            +cff3*psi3(j,i)+cff4*psi4(j,i)
         enddo
      enddo

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine phi_from_psi(ii,jj,psi1,psi2,phi1,phi2)
      implicit none
      integer ii,jj,i,j,leapfrog
      real*8 psi1(ii,jj,2),psi2(ii,jj,2),phi1(ii,jj,2),phi2(ii,jj,2)
     & ,theta11,theta12,theta21,theta22

      common /FOUR/ theta11,theta12,theta21,theta22

      do leapfrog=1,2
         do j=1,jj
            do i=1,ii
               phi1(i,j,leapfrog)= theta11*psi1(i,j,leapfrog)
     &                            +theta12*psi2(i,j,leapfrog)
               phi2(i,j,leapfrog)= psi1(i,j,leapfrog)
     &                            -psi2(i,j,leapfrog)
            enddo
         enddo
      enddo

      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine psi_from_phi(ii,jj,phi1,phi2,psi1,psi2)
      implicit none
      integer ii,jj,i,j,leapfrog
      real*8 psi1(ii,jj,2),psi2(ii,jj,2)
     & ,phi1(ii,jj,2),phi2(ii,jj,2)
     & ,omega11,omega12,omega21,omega22

      common /FIVE/ omega11,omega12,omega21,omega22

      do leapfrog=1,2
         do j=1,jj
            do i=1,ii
               psi1(i,j,leapfrog)=         phi1(i,j,leapfrog)
     &                            +omega12*phi2(i,j,leapfrog)
               psi2(i,j,leapfrog)=         phi1(i,j,leapfrog)
     &                            +omega22*phi2(i,j,leapfrog)
            enddo
         enddo
      enddo

      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine z_from_phi_and_rel(ii,jj,SS,phi1,phi2,rel1,rel2,z1,z2)
      implicit none
      integer ii,jj,i,j,leapfrog
      real*8 phi1(ii,jj,2),phi2(ii,jj,2),rel1(ii,jj,2),rel2(ii,jj,2)
     & ,z1(ii,jj,2),z2(ii,jj,2)
     & ,SS

      do leapfrog=1,2
         do j=1,jj
            do i=1,ii
               z1(i,j,leapfrog)=rel1(i,j,leapfrog)
               z2(i,j,leapfrog)=rel2(i,j,leapfrog)
     &                                     -SS*phi2(i,j,leapfrog)
            enddo
         enddo
      enddo

      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine rel_from_psi(islip,ii,jj,phi1,phi2,rel1,rel2)
      implicit none
      integer jj,ii,islip,leapfrog,j,i
      real*8 phi1(ii,jj,2),phi2(ii,jj,2),rel1(ii,jj,2),rel2(ii,jj,2)

      do leapfrog=1,2
         do j=2,jj-1
            do i=2,ii-1
               rel1(i,j,leapfrog)=-4.*phi1(i,j,leapfrog)
     &                        +phi1(i-1,j,leapfrog)+phi1(i,j+1,leapfrog)
     &                        +phi1(i+1,j,leapfrog)+phi1(i,j-1,leapfrog)
               rel2(i,j,leapfrog)=-4.*phi2(i,j,leapfrog)
     &                        +phi2(i-1,j,leapfrog)+phi2(i,j+1,leapfrog)
     &                        +phi2(i+1,j,leapfrog)+phi2(i,j-1,leapfrog)
            enddo
         enddo

         i=1
         do j=2,jj-1
            rel1(i,j,leapfrog)=-4.*phi1(i,j,leapfrog)
     &                        +phi1( ii,j,leapfrog)+phi1(i,j+1,leapfrog)
     &                        +phi1(i+1,j,leapfrog)+phi1(i,j-1,leapfrog)
            rel2(i,j,leapfrog)=-4.*phi2(i,j,leapfrog)
     &                        +phi2( ii,j,leapfrog)+phi2(i,j+1,leapfrog)
     &                        +phi2(i+1,j,leapfrog)+phi2(i,j-1,leapfrog)
         enddo

         i=ii
         do j=2,jj-1
            rel1(i,j,leapfrog)=-4.*phi1(i,j,leapfrog)
     &                        +phi1(i-1,j,leapfrog)+phi1(i,j+1,leapfrog)
     &                        +phi1(  1,j,leapfrog)+phi1(i,j-1,leapfrog)
            rel2(i,j,leapfrog)=-4.*phi2(i,j,leapfrog)
     &                        +phi2(i-1,j,leapfrog)+phi2(i,j+1,leapfrog)
     &                        +phi2(  1,j,leapfrog)+phi2(i,j-1,leapfrog)
         enddo
         call boundary_condition(islip,leapfrog,ii,jj
     &                                ,rel1,rel2,phi1,phi2)
      enddo

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine zeta_from_psi_and_rel(ii,jj,S1,S2,psi1,psi2,rel1,rel2
     &                                                    ,zeta1,zeta2)
      implicit none
      integer ii,jj,i,j,leapfrog
      real*8 psi1(ii,jj,2),psi2(ii,jj,2),rel1(ii,jj,2),rel2(ii,jj,2)
     & ,zeta1(ii,jj,2),zeta2(ii,jj,2)
     & ,S1,S2

      do leapfrog=1,2
      do j=1,jj
         do i=1,ii
            zeta1(i,j,leapfrog)=rel1(i,j,leapfrog)
     &                 -S1*(psi1(i,j,leapfrog)-psi2(i,j,leapfrog))
            zeta2(i,j,leapfrog)=rel2(i,j,leapfrog)
     &                 -S2*(psi2(i,j,leapfrog)-psi1(i,j,leapfrog))
         enddo
      enddo
      enddo

      return
      end
c
c----------------------------------------------------------------------
c
      subroutine zeta_interior(leapfrog,m,T,ii,jj,S1,S2,U_0
     &                                ,beta1_y,beta2_y
     &                                ,bot
     &                                ,visc,visc_bot
     &                                ,psi1,psi2,rel1,rel2
     &                                ,pv1,pv2,pv1_x,q1,q2,u1,u2,v1,v2
     &                                ,diss1,diss2
     &                                ,zeta1,zeta2
     &                                ,D1_s,D2_s,D1_n,D2_n
     &                                ,A1_s,A2_s,A1_n,A2_n)
      implicit none
      integer ii,jj,i,j,leapfrog,l,m
     & ,i1,j1,im1,ip1,jm1,jp1,im2,ip2,jm2,jp2
      real*8 psi1(ii,jj,2),psi2(ii,jj,2)
     & ,beta1_y(jj),beta2_y(jj),bot(0:ii+1,jj)
     & ,rel1(ii,jj,2),rel2(ii,jj,2),zeta1(ii,jj,2),zeta2(ii,jj,2)

     & ,pv1(0:ii+1,jj),pv2(0:ii+1,jj),pv1_x(ii,jj)
     & ,q1(0:ii+1,jj),q2(0:ii+1,jj)
     & ,u1(0:ii+1,jj),u2(0:ii+1,jj)
     & ,v1(0:ii+1,jj),v2(0:ii+1,jj)
     & ,diss1(ii,jj),diss2(ii,jj)

     & ,T,S1,S2,U_0,visc,visc_bot
     & ,cff1,cff2
     & ,theta11,theta12,theta21,theta22
     & ,rho1_11,rho1_12,rho1_21,rho1_22
     & ,rho2_11,rho2_12,rho2_21,rho2_22
     & ,D1_s,D2_s,D1_n,D2_n
     & ,A1_s,A2_s,A1_n,A2_n

      common /FOUR/ theta11,theta12,theta21,theta22
      common /SIX/ rho1_11,rho1_12,rho1_21,rho1_22
     &            ,rho2_11,rho2_12,rho2_21,rho2_22

c
c--- auxiliary arrays
c
      do j=1,jj
         do i=1,ii
            pv1(i,j)=rel1(i,j,m)-S1*(psi1(i,j,m)-psi2(i,j,m))
            pv2(i,j)=rel2(i,j,m)-S2*(psi2(i,j,m)-psi1(i,j,m))
         enddo
         pv1(0,j)=pv1(ii,j)
         pv2(0,j)=pv2(ii,j)
         pv1(ii+1,j)=pv1(1,j)
         pv2(ii+1,j)=pv2(1,j)
      enddo

      do j=1,jj
         do i=0,ii+1
            q1(i,j)=pv1(i,j)+beta1_y(j)
            q2(i,j)=pv2(i,j)+beta2_y(j)+bot(i,j)
         enddo
      enddo

      do j=1,jj
         do i=2,ii-1
            pv1_x(i,j)=.5*(pv1(i+1,j)-pv1(i-1,j))
         enddo
         pv1_x(1,j)=.5*(pv1(2,j)-pv1(ii,j))
         pv1_x(ii,j)=.5*(pv1(1,j)-pv1(ii-1,j))
      enddo

      do j=2,jj-1   ! 2*velocity
         do i=1,ii
            u1(i,j)=psi1(i,j-1,m)-psi1(i,j+1,m)
            u2(i,j)=psi2(i,j-1,m)-psi2(i,j+1,m)
         enddo
         u1(0,j)=u1(ii,j)
         u2(0,j)=u2(ii,j)
         u1(ii+1,j)=u1(1,j)
         u2(ii+1,j)=u2(1,j)
      enddo

      do j=2,jj-1   ! 2*velocity
         do i=2,ii-1
            v1(i,j)=psi1(i+1,j,m)-psi1(i-1,j,m)
            v2(i,j)=psi2(i+1,j,m)-psi2(i-1,j,m)
         enddo
         v1(1,j)=psi1(2,j,m)-psi1(ii,j,m)
         v2(1,j)=psi2(2,j,m)-psi2(ii,j,m)
         v1(ii,j)=psi1(1,j,m)-psi1(ii-1,j,m)
         v2(ii,j)=psi2(1,j,m)-psi2(ii-1,j,m)
         v1(0,j)=v1(ii,j)
         v2(0,j)=v2(ii,j)
         v1(ii+1,j)=v1(1,j)
         v2(ii+1,j)=v2(1,j)
      enddo

      do j=2,jj-1
         do i=2,ii-1
            diss1(i,j)=visc*(-4.*rel1(i,j,leapfrog)
     &                  +rel1(i-1,j,leapfrog)+rel1(i+1,j,leapfrog)
     &                  +rel1(i,j-1,leapfrog)+rel1(i,j+1,leapfrog)
     &                      )
            diss2(i,j)=visc*(-4.*rel2(i,j,leapfrog)
     &                  +rel2(i-1,j,leapfrog)+rel2(i+1,j,leapfrog)
     &                  +rel2(i,j-1,leapfrog)+rel2(i,j+1,leapfrog)
     &                      )
         enddo
         diss1(1,j)=visc*(-4.*rel1(1,j,leapfrog)
     &                  +rel1(ii, j,leapfrog)+rel1(2,j,leapfrog)
     &                  +rel1(1,j-1,leapfrog)+rel1(1,j+1,leapfrog)
     &                   )
         diss2(1,j)=visc*(-4.*rel2(1,j,leapfrog)
     &                  +rel2(ii, j,leapfrog)+rel2(2,j,leapfrog)
     &                  +rel2(1,j-1,leapfrog)+rel2(1,j+1,leapfrog)
     &                   )
         diss1(ii,j)=visc*(-4.*rel1(ii,j,leapfrog)
     &                  +rel1(ii-1,j,leapfrog)+rel1(1,j,leapfrog)
     &                  +rel1(ii,j-1,leapfrog)+rel1(ii,j+1,leapfrog)
     &                    )
         diss2(ii,j)=visc*(-4.*rel2(ii,j,leapfrog)
     &                  +rel2(ii-1,j,leapfrog)+rel2(1,j,leapfrog)
     &                  +rel2(ii,j-1,leapfrog)+rel2(ii,j+1,leapfrog)
     &                    )
      enddo
c
c--- L1
c
      do j=2,jj-1
         do i=1,ii
            zeta1(i,j,leapfrog)=zeta1(i,j,leapfrog)+T*(
     & -.25*(
     &       u1(i+1,j)*q1(i+1,j)-u1(i-1,j)*q1(i-1,j)
     &      +v1(i,j+1)*q1(i,j+1)-v1(i,j-1)*q1(i,j-1)
     &      )
     & -U_0*pv1_x(i,j)
     & +diss1(i,j)
     &                                                )
         enddo
      enddo
c
c--- L2
c
      do j=2,jj-1
         do i=1,ii
            zeta2(i,j,leapfrog)=zeta2(i,j,leapfrog)+T*(
     & -.25*(
     &       u2(i+1,j)*q2(i+1,j)-u2(i-1,j)*q2(i-1,j)
     &      +v2(i,j+1)*q2(i,j+1)-v2(i,j-1)*q2(i,j-1)
     &      )
     & +diss2(i,j)
     & -visc_bot*rel2(i,j,leapfrog)
     &                                                )
         enddo
      enddo

      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine mean_periodic(ii,jj,psi,psiC)
      implicit none
      integer ii,jj,i,j
      real*8 psi(ii,jj),psiC,sum

      sum=0.
      do j=2,jj-1
         do i=1,ii
            sum=sum+psi(i,j)
         enddo
      enddo
      do i=1,ii
         sum=sum+0.5*(psi(i,1)+psi(i,jj))
      enddo
      psiC=sum/dfloat(ii*(jj-1))

      return
      end
c
c-------------------------------------------------------------------------
c
      subroutine subtract_mean(ii,jj,psi)
      implicit none      
      integer ii,jj,ii1,jj1,i,j
      real psi(ii,jj)
      real*8 psiC,sum

      ii1=ii-1
      jj1=jj-1
      sum=0.
      do j=2,jj1
         do i=2,ii1
            sum=sum+psi(i,j)
         enddo                 
      enddo           
      do j=2,jj1
         sum=sum+0.5*(psi(1,j)+psi(ii,j))
      enddo
      do i=2,ii1
         sum=sum+0.5*(psi(i,1)+psi(i,jj))
      enddo                    
      sum=sum+psi(1,1)  ! 0.25*(* + * + * + *)
      psiC=sum/float(ii1*jj1)

      do j=1,jj
         do i=1,ii
            psi(i,j)=psi(i,j)-psiC
         enddo
      enddo                    

      return
      end
c
c--------------------------------------------------------------------------
c
      subroutine subtract_mean_periodic(ii,jj,psi)
      implicit none
      integer ii,jj,i,j
      real*8 psi(ii,jj),psiC,sum

      sum=0.
      do j=2,jj-1
         do i=1,ii
            sum=sum+psi(i,j)
         enddo
      enddo
      do i=1,ii
         sum=sum+0.5*(psi(i,1)+psi(i,jj))
      enddo
      psiC=sum/dfloat(ii*(jj-1))

      do j=1,jj
         do i=1,ii
            psi(i,j)=psi(i,j)-psiC
         enddo
      enddo

      return
      end
c
c--------------------------------------------------------------------------
c
      subroutine energy(H1,H2,S1,S2,ii,jj,psi1,psi2,ekin1,ekin2,epot)
      implicit none
      integer ii,jj,i,j,ii1,jj1
      real*8 H1,H2,S1,S2,ekin1,ekin2,epot,H
     & ,psi1(ii,jj,2),psi2(ii,jj,2)
     & ,fac1,fac2,cff

      ii1=ii-1
      jj1=jj-1
      H=H1+H2
      ekin1=0.
      ekin2=0.
      epot=0.

      do j=2,jj1
         do i=2,ii1
            epot=epot+(psi1(i,j,2)-psi2(i,j,2))**2
         enddo
      enddo
      epot=epot+float(ii1+jj1-1)*(psi1(1,1,2)-psi2(1,1,2))**2
      do j=2,jj
         do i=2,ii
           ekin1=ekin1+.25*(
     &   (psi1(i-1,j,2)+psi1(i,j,2)-psi1(i-1,j-1,2)-psi1(i,j-1,2))**2
     &  +(psi1(i,j-1,2)+psi1(i,j,2)-psi1(i-1,j-1,2)-psi1(i-1,j,2))**2
     &                     )
           ekin2=ekin2+.25*(
     &   (psi2(i-1,j,2)+psi2(i,j,2)-psi2(i-1,j-1,2)-psi2(i,j-1,2))**2
     &  +(psi2(i,j-1,2)+psi2(i,j,2)-psi2(i-1,j-1,2)-psi2(i-1,j,2))**2
     &                     )
         enddo
      enddo
      cff=0.5/dfloat(ii1*jj1)
      ekin1=cff*(H1/H)*ekin1
      ekin2=cff*(H2/H)*ekin2
      epot =0.5*cff*((S1*H1+S2*H2)/H)*epot

      return
      end
c
c-------------------------------------------------------------------------
c
      function ran1(idum)
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836
     & ,ntab=32,ndiv=1+(im-1)/ntab,eps=1.2E-7,rnmx=1.-eps)
      integer iv(ntab)
      save iv,iy
      data iv /ntab*0/, iy /0/

      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          if (idum.lt.0) idum=idum+im
          if (j.le.ntab) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)

      return
      end

#ifdef advect_particles
c
c-----------------------------------------------------------------------------
c
      subroutine deriv(x_c,y_c,psi_c,ii,jj,x,y,u,v)
      implicit none
      integer ii,jj,i,j
      real*8 dx,cff
      parameter(dx=.05D0,cff=.5D0/dx)
      real*8 psi_c(0:ii+2,jj),x_c(0:ii+2),y_c(jj)
     & ,x,y,u,v,x1,x2,y1,y2,p1,p2,p3,p4

      x1=x-dx
      x2=x+dx
      y1=y-dx
      y2=y+dx
      call cubic(x_c,y_c,psi_c,ii,jj,x1,y,p1)
      call cubic(x_c,y_c,psi_c,ii,jj,x,y2,p2)
      call cubic(x_c,y_c,psi_c,ii,jj,x2,y,p3)
      call cubic(x_c,y_c,psi_c,ii,jj,x,y1,p4)
      u=cff*(p4-p2)
      v=cff*(p3-p1)
               
      return
      end
c
c-----------------------------------------------------------------------------
c
      subroutine cubic(x_c,y_c,f_c,ii,jj,x,y,f)
      implicit none
      integer ii,jj,i,j,ic,jc
      real*8 c1,c2,dx,dy
      parameter(c1=1.D0/6.D0,c2=2.D0/6.D0)
      real*8 f_c(0:ii+2,jj),x_c(0:ii+2),y_c(jj),cff(4),f_aux(4)
     & ,x,y,f,ax,ax2,ax3,ay,ay2,ay3

      dx=x_c(5)-x_c(4)
      dy=y_c(5)-y_c(4)

      ic=int(x)+1
      do i=1,ii-1
         if (x.ge.x_c(i) .and. x.lt.x_c(i+1)) then
           ic=i
         endif
      enddo
      if(ic.ge.ii-1) ic=ii-2

      jc=int(y)+1
      do j=1,jj-1
         if (y.ge.y_c(j) .and. y.lt.y_c(j+1)) then
           jc=j
         endif 
      enddo 
      if(jc.le.1) jc=2
      if(jc.ge.jj-1) jc=jj-2

      ax=(x-x_c(ic))/dx
      ay=(y-y_c(jc))/dy
      ax2=ax**2
      ax3=ax**3
      ay2=ay**2
      ay3=ay**3
c
c--- y-interpolation
c
      cff(1)=  -c2*ay+.5*ay2-c1*ay3
      cff(2)=1.-.5*ay   -ay2+.5*ay3
      cff(3)=      ay+.5*ay2-.5*ay3
      cff(4)=  -c1*ay       +c1*ay3
      do j=1,4
         f_aux(j)=0.
         do i=1,4
            f_aux(j)=f_aux(j)+cff(i)*f_c(ic-2+i,jc-2+j)
         enddo
      enddo
c
c--- x-interpolation
c
      cff(1)=  -c2*ax+.5*ax2-c1*ax3
      cff(2)=1.-.5*ax   -ax2+.5*ax3
      cff(3)=      ax+.5*ax2-.5*ax3
      cff(4)=  -c1*ax       +c1*ax3
      f= cff(1)*f_aux(1)+cff(2)*f_aux(2)
     &  +cff(3)*f_aux(3)+cff(4)*f_aux(4)

      return
      end


c
c-----------------------------------------------------------------------------
c
      subroutine coarsen(f_f,f_c,ii,jj,iic,jjc,ncoarse)
c
      integer jfi, ifi
      real*8 f_c(0:iic+2,jjc), f_f(0:ii+2,jj)
c
      do j=1,jjc
         jfi = ncoarse*(j-1) + 1 
         do i=1,iic
            ifi = ncoarse*(i-1) + 1
            f_c(i,j) = f_f(ifi,jfi)
         enddo
         f_c(0,j) = f_c(iic,j)
         f_c(iic+1,j) = f_c(1,j)
         f_c(iic+2,j) = f_c(2,j)
      enddo
c
      return
      end
c
#endif
c---------------------------------------------------------------------
