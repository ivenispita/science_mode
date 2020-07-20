function [Vo,En,U,V,ETA] = isw_2d_cb(dx,dy,Lx,Ly,H,g,dt,ts,u0,v0,eta0)
	%
	% Integrate the 2D-inviscid shallow water model equations,
	% for a closed rectangular basin spatially discretized following 
	% a C-grid scheme. 
	%
	% This simple model is integrated in the dimensional form for a non-rotating frame !
	% 
	% The models is unforced and require initial conditions
	%
	%
	% Input:
	%	dx, dy: spatial resolution
	%
	%	Lx, Ly: width and height of the basin 
	%
	%	H: topography (i. e., local depth)
	%
	%	g: acceleration due to gravity
	%
	%	dt: time step
	%
	%	ts: Saving outputs time-steps
	%
	%	u0, v0, eta0: initial condition
	%
	%
	% Returns:
	%	Vo: Total volume of displaced water
	%
	%	En: Total energy of the system
	%
	%	U,V,ETA: Solutions at ts
	%
	%
	% Author: Tiago Bilo
	% CFD - Fall 2016
	% Problem set 4:
	% 2. Numerical Experimentation	 


    % Grid set-up
    [x,xe,Nx] = FDGrid(0,Lx,dx); 
    [y,ye,Ny] = FDGrid(0,Ly,dy);

    [x_eta,y_eta] = meshgrid(x,y);      % eta,H grid points
    [x_u,y_u] = meshgrid(xe,y);         % u,U   grid points
    [x_v,y_v] = meshgrid(x,ye);         % v,V   grid points
    [xe,ye] = meshgrid(xe,ye);          % vorticiy grid points (or cell edges)

    size_ts = length(ts);
    size_u = size(x_u);
    size_v = size(x_v);
    size_eta = size(x_eta);


    % Initializing the variables
    Vo = ones(1,size_ts)*nan;
    En = ones(1,size_ts)*nan;
    U = ones(size_u(1),size_u(2),size_ts)*nan;
    V = ones(size_v(1),size_v(2),size_ts)*nan;
    ETA = ones(size_eta(1),size_eta(2),size_ts)*nan;

    % Initial conditions
    u = u0;
    v = v0;
    eta = eta0;

    % Integration
    for t = [dt:dt:ts(end)]

    	% Updating the variables
		[u,v,eta] = rk3_isw(u,v,eta,dx,dy,dt,H,g);			% RK3 time-stepping scheme

		% Enforcing boundary conditions
		% Closed boundaries: No normal flow at the boundaries
		u(:,[1 end]) = 0.0;
		v([1 end],:) = 0.0;


		% Total volume of displaced water
		volume = sum(sum(eta))*dx*dy;

		% Energy of the system
		% x and y axis corresponds to the 2nd and 1st dimensions of the matrices, respectively !
		u_xaverage = yop2_2d(eta,u,0.0,[0.5 0.5]);
		v_yaverage = xop2_2d(eta,v,0.0,[0.5 0.5]);

		energy = (g.*eta.*eta)+(H.*u_xaverage.*u_xaverage) ...
			+(H.*v_yaverage.*v_yaverage);

		energy = sum(sum(energy))*dx*dy/2.0;

		i_ts = find(ts == t);
		if length(i_ts) == 1
			Vo(i_ts) = volume;
			En(i_ts) = energy;
			U(:,:,i_ts) = u;
			V(:,:,i_ts) = v;
			ETA(:,:,i_ts) = eta;
		end 
    end



