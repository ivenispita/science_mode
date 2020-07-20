function [Vo,En,Et,U,V,ETA,ZETA] = nl_isw_2d_cb(dx,dy,xlim,ylim,f,a,H,g,dt,ts,u0,v0,eta0,mask_u,mask_v)
	%
	% Integrate the 2D-inviscid non-linear shallow water model equations,
	% for a closed rectangular basin spatially discretized following 
	% an Arakawa C-grid scheme. 
	%
	%
	% The models is unforced and require initial conditions
	%
	%
	% Input:
	%	dx, dy: spatial resolution
	%
	%	xlim,ylim: physical limits of the basin
	%
	%	f: planetary vorticity  
	%
	%	H: topography (i. e., local depth)
	%
	%	a: non-linear terms parameter. If a = 1 - Include the non-linear terms;
	%	if a = 0 - Solve the linear system
	%
	%	g: acceleration due to gravity
	%
	%	dt: time step
	%
	%	ts: Saving outputs time-steps
	%
	%	u0, v0, eta0: initial condition
	%
	%   mask_u, mask_v: Topographic masks
	%
	% Returns:
	%	Vo: Total volume of displaced water
	%
	%	En: Total energy of the system
	%
	%	Et: Potential enstrophy of the system
	%
	%	U,V,ETA,ZETA: Solutions at ts
	%
	%
	% Author: Tiago Bilo
	% CFD - Fall 2016
	% Class Project
 


    % Grid set-up
    [x,xe,Nx] = FDGrid(xlim(1),xlim(2),dx); 
    [y,ye,Ny] = FDGrid(ylim(1),ylim(2),dy);

    [x_eta,y_eta] = meshgrid(x,y);      % eta,H grid points
    [x_u,y_u] = meshgrid(xe,y);         % u,U   grid points
    [x_v,y_v] = meshgrid(x,ye);         % v,V   grid points
    [x_q,y_q] = meshgrid(xe,ye);        % vorticiy grid points (or cell edges)

    size_ts = length(ts);
    size_u = size(x_u);
    size_v = size(x_v);
    size_eta = size(x_eta);
	size_q = size(x_q);

    % Initializing the variables
    Vo = ones(1,size_ts)*nan;
    En = ones(1,size_ts)*nan;
    Et = ones(1,size_ts)*nan;    
    U = ones(size_u(1),size_u(2),size_ts)*nan;
    V = ones(size_v(1),size_v(2),size_ts)*nan;
    ZETA = ones(size_q(1),size_q(2),size_ts)*nan;    
    ETA = ones(size_eta(1),size_eta(2),size_ts)*nan;

    % Initial conditions
    u = u0;
    v = v0;
    eta = eta0;

    % Integration
    for t = [0:1:ts(end)]

    	% Updating the variables
		[u,v,eta] = rk3_isw(u,v,eta,dx,dy,dt,f,a,H,g,mask_u,mask_v);			% RK3 time-stepping scheme

		% Enforcing boundary conditions
		% Closed boundaries: No normal flow at the boundaries
		u(:,[1 end]) = 0.0;
		v([1 end],:) = 0.0;

		u = u.*mask_u;
		v = v.*mask_v;

		% Total volume of displaced water
		volume = sum(sum(eta))*dx*dy;

		%% Energy of the system
		% x and y axis corresponds to the 2nd and 1st dimensions of the matrices, respectively !
		h = H+eta;												% Total water column height
		u2 = u.*u;
		v2 = v.*v;

		u2_xaverage = yop2_2d(eta,u2,0.0,[0.5 0.5]);
		v2_yaverage = xop2_2d(eta,v2,0.0,[0.5 0.5]);

		energy = (g.*eta.*eta)+h.*(u2_xaverage+v2_yaverage);
		energy = sum(sum(energy))*dx*dy/2.0;


		%% Enstrophy of the system
		% x and y axis corresponds to the 2nd and 1st dimensions of the matrices, respectively !
		h_average = zeros(size_q);
		ha = yop2_2d(u,h,0.0,[0.5 0.5]); 
		h_average(:,2:end) = xop2_2d(v,ha(:,2:end),0.0,[0.5 0.5]);

		% Enforcing that potential vorticity is not defined at the boundaries
		% i. e., f/h is not defined
		h_average([1 end],:) = nan;
		h_average(:,[1 end]) = nan;

		dvdx = zeros(size_q);
		dudy = zeros(size_q);

		dvdx(2:end,:) = yop2_2d(u,v(2:end,:),0.0,[1/dx -1/dx]);
		dudy(:,2:end) = xop2_2d(v,u(:,2:end),0.0,[1/dy -1/dy]);

		zeta = dvdx-dudy;
		q = (zeta+f)./h_average;

		enstrophy = q.*q;
		enstrophy = nansum(nansum(enstrophy))*dx*dy/2.0;

		i_ts = find(ts == t);
		if length(i_ts) == 1
			Vo(i_ts) = volume;
			En(i_ts) = energy;
			Et(i_ts) = enstrophy;			
			U(:,:,i_ts) = u;
			V(:,:,i_ts) = v;
			ZETA(:,:,i_ts) = zeta;
			ETA(:,:,i_ts) = eta;
		end 
    end



