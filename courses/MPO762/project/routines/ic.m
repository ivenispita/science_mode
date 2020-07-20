function [u,v,eta] = ic(x,xe,y,ye,f,g,x0,y0,U,Rd)
	%
	% Create anticyclonic eddies
	% The eddies are in geostrophic balance 
	% 
	% xe,ye: Grid cell edges coordinates
	% x,y: Grid cell center coordinates
	% x0, y0: Core coordinates
	% f: Coriolis Parameter
	% g: acceleration due to gravity or reduced gravity term
	% Rd: velocity decay scale 
	% U: Maximum velocity scale
	%
	%
	% Author: Tiago Bilo
	% CFD - Fall 2016
	% Class project

	% Grid points
	[x_eta,y_eta] = meshgrid(x,y);      % eta   grid points
	[x_u,y_u] = meshgrid(xe,y);         % u     grid points
	[x_v,y_v] = meshgrid(x,ye);         % v     grid points

	% Gaussian structure
	phi_v = -((x_v-x0).^2.0)/(Rd.^2);
	gamma_v = -((y_v-y0).^2.0)/(Rd.^2);

	gamma_u = -((y_u-y0).^2.0)/(Rd.^2);
	phi_u = -((x_u-x0).^2.0)/(Rd.^2);

	gamma_eta = -((y_eta-y0).^2.0)/(Rd.^2);
	phi_eta = -((x_eta-x0).^2.0)/(Rd.^2);


	% Velocity fields
	if f>0
		vf = -1.0;
		uf = 1.0;
		etaf = 1.0;
	else
		vf = 1.0;
		uf = -1.0;
		etaf = -1.0;		
	end		


	XYv = x_v-x0;
	v = vf*(U.*XYv.*exp(phi_v).*exp(gamma_v)./Rd);

	XYu = y_u-y0;
	u = uf*(U.*XYu.*exp(phi_u).*exp(gamma_u)./Rd);

	eta = etaf*(U.*f.*Rd.*exp(phi_eta).*exp(gamma_eta)/(2*g));
