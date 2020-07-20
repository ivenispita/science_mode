function [mask_u,mask_v,mask_q] = build_ridge(x,xe,x0,y,ye,L,d,r)
	% 
	% Build a meridional ridge based on its radius r    
	% and distance between them d at the longitude  
	% x0
	%
	% Author: Tiago Bilo
	% CFD - Fall 2016
	% Class project

	y0 = [-L/2.0:(2*r)+d:L/2.0];

	% Grid points
	[x_u,y_u] = meshgrid(xe,y);         % u     grid points
	[x_v,y_v] = meshgrid(x,ye);         % v     grid points
	[x_q,y_q] = meshgrid(xe,ye);        % q     grid points

	% Masks 
	mask_u = ones(size(x_u));
	mask_v = ones(size(x_v));
	mask_q = ones(size(x_q));


	for i = [1:length(y0)]
		rru = ((x_u-x0).^2) + ((y_u-y0(i)).^2);
		rrv = ((x_v-x0).^2) + ((y_v-y0(i)).^2);
		rrq = ((x_q-x0).^2) + ((y_q-y0(i)).^2);

		iu = rru<=r.^2;
		iv = rrv<=r.^2;
		iq = rrq<=r.^2;

		mask_u(iu) = 0.0;
		mask_v(iv) = 0.0;
		mask_q(iq) = 0.0;		
	end

