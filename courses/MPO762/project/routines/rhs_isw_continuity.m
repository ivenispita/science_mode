function r = rhs_isw_continuity(u,v,eta,dx,dy,H)
	%
	% Compute time-tendency (i. e., the convergence of the mass flux) 
	% of the sea-surface height
	%
	% usage is
	% r = rhs_isw_continuity(u,v,eta,dx,dy,H)
	% where u,v: are the velocity components
	%
	%       eta: is the sea-surface height field at the current time step
	% 
	%       dx,dy: are the grid spacing
	%
	%       H: is the topography (i. e., local depth)
	%
	% Author: Tiago Bilo
	% CFD - Fall 2016
	% Class Project

	% Verifying size of H (if H = cte, it will be expanded to all grid cell centers)
	if length(H) == 1
		H = ones(size(eta))*H;
	end

	% Total water column height
	h = H+eta;

	% Taking the average depth at u, v grid points
	% x and y axis corresponds to the 2nd and 1st dimensions of the matrices, respectively !
	hu = yop2_2d(u,h,0.0,[0.5 0.5]);  
	hv = xop2_2d(v,h,0.0,[0.5 0.5]);		

	% Mass fluxes U,V far from the boundaries (at boundaries = 0.0)
	U = u.*hu;
	V = v.*hv;

	% Differentiate U and V
	dU = yop2_2d(eta,U,0.0,[1.0/dx -1.0/dx]);
	dV = xop2_2d(eta,V,0.0,[1.0/dy -1.0/dy]);

	% Right Hand Side of the continuity: i. e., mass convergence term
	r = -dU-dV;
