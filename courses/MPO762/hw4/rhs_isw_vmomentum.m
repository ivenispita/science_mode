function r = rhs_isw_vmomentum(v,eta,dy,g)
	%
	% Compute time-tendency (i. e., meridional pressure gradient) 
	% of the meridional velocity
	%
	% usage is
	% r = rhs_isw_umomentum(v,eta,dy,g)
	% where v: is the meridional velocity component
	%
	%       eta: is the sea-surface height field at the current time step
	% 
	%       dy: is the meridional grid spacing
	%
	%		g: acceleration due to gravity
	%
	% Author: Tiago Bilo
	% CFD - Fall 2016
	% Problem set 4:
	% 2. Numerical Experimentation


	dndy = xop2_2d(v,eta,0.0,[1.0/dy -1.0/dy]);
	r = -g*dndy;

	% Enforcing boundary conditions
	% Closed boundaries: No normal flow at the boundaries (i. e., no v tendencies)
	r([1 end],:) = 0.0;