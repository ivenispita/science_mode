function [uu,vu,etau] = rk3_isw(u,v,eta,dx,dy,dt,f,a,H,g,mask_u,mask_v)
	%
	% RK3 Scheme to time steps the variables u, v and eta through a single step of size depth
	%
	% usage is
	% [uu,vu,etau] = rk3_isw(u,v,eta,dx,dy,dt,f,a,H,g,mask_u,mask_v)
	%
	% where u,v,eta: are the variables of the inviscid shallow water model 
	%
	%       dx,dy: are the grid spacing
	%
	%       dt: is the time step
	%
	%		f: is the planetary vorticity
	%
	%		a: non-linear terms parameter. If a = 1 - Include the non-linear terms;
	%		if a = 0 - Solve the linear system
	%
	%       H: is the topography (i. e., local depth)
	%
	%       g: acceleration due to gravity
	%
	%       mask_u, mask_v: Topographic masks
	%
	% Author: Tiago Bilo
	% CFD - Fall 2016
	% Class Project




	% Eta update (i. e., integrate the continuity equation)
	% Stage 1 update
	r = rhs_isw_continuity(u,v,eta,dx,dy,H);
	etat1 = eta + dt*r;               			

	r = rhs_isw_umomentum(u,v,eta,dx,dy,a,f,g,H,mask_u);
	ut1 = u + dt*r;               			    


	r = rhs_isw_vmomentum(u,v,eta,dx,dy,a,f,g,H,mask_v);
	vt1 = v + dt*r;               			    


	% Stage 2 update
	r = rhs_isw_continuity(ut1,vt1,etat1,dx,dy,H);
	etat2 = 0.75*eta + 0.25*(etat1 + dt*r); 		

	r = rhs_isw_umomentum(ut1,vt1,etat1,dx,dy,a,f,g,H,mask_u);
	ut2 = 0.75*u + 0.25*(ut1 + dt*r); 		    

	r = rhs_isw_vmomentum(ut1,vt1,etat1,dx,dy,a,f,g,H,mask_v);
	vt2 = 0.75*v + 0.25*(vt1 + dt*r); 		    


	% Stage 3 final update
	r = rhs_isw_continuity(ut2,vt2,etat2,dx,dy,H);
	etau = (eta + 2.0*(etat2 + dt*r))/3.0;  		

	r = rhs_isw_umomentum(ut2,vt2,etat2,dx,dy,a,f,g,H,mask_u);
	uu = (u + 2.0*(ut2 + dt*r))/3.0;  		    

	r = rhs_isw_vmomentum(ut2,vt2,etat2,dx,dy,a,f,g,H,mask_v);
	vu = (v + 2.0*(vt2 + dt*r))/3.0;  		    

