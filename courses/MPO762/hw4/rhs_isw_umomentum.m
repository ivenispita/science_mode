function r = rhs_isw_umomentum(u,eta,dx,g)
	%
	% Compute time-tendency (i. e., zonal pressure gradient) 
	% of the zonal velocity
	%
	% usage is
	% r = rhs_isw_umomentum(u,eta,dx,g)
	% where u: is the zonal velocity component
	%
	%       eta: is the sea-surface height field at the current time step
	% 
	%       dx: is the zonal grid spacing
	%
	%		g: acceleration due to gravity
	%
	% Author: Tiago Bilo
	% CFD - Fall 2016
	% Problem set 4:
	% 2. Numerical Experimentation


	dndx = yop2_2d(u,eta,0.0,[1.0/dx -1.0/dx]);
	r = -g*dndx;

	% Enforcing boundary conditions
	% Closed boundaries: No normal flow at the boundaries (i. e., no u tendencies)
	r(:,[1 end]) = 0.0;
