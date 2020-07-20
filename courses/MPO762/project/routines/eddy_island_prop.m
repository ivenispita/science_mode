% Author: Tiago Bilo
% CFD - Fall 2016
% Class project:
% Self-propagating eddies and obstacles interactions in a zonal channel 
% solved by the Inviscid Shallow Water equations 

clear all 
close all 

disp(['Model options'])
disp(' ')
h = input('Hemisphere: HS = -1 and HN = 1  ' );
exper = input('Topography [y/n]: ');

if h < 0 
	hm = 'hs';
else
	hm = 'hn';
end


%%%% Model parameters 
% Physical parameters
H = 1000.0;                       % Cte depth (1000 m)
g = 0.0045;                       % Reduced gravity gravity (<< 10 m/s2 and < 0) or acceleration due to gravity (~10 m/s2)
f0 = gsw_f(h*13);		  		  % Coriolis parameter (s-1)
Rd = sqrt(H*g)/abs(f0);			  % Rossby Radius of deformation

Lx = 1500000.0;				      % Channel length (m) 
Ly = 600000.0;                    % Channel width (m)


% I have to figure out the dimensions of the islands
% and the distance between them 

% Grid parameters (Arakawa C-grid)
dx = 2000;            		  			        % zonal resolution 2 km (m)
dy = dx;                          				% meridional resolution (m)
dt = 500;                       	            % Time step of ~8 min (s) 

xlim = [-Lx/2.0,Lx/2.0]; 		  				% Channel zonal boundaries              
ylim = [-Ly/2.0,Ly/2.0];		  				% Channel meridional boundaries


[x,xe,Nx] = FDGrid(xlim(1),xlim(2),dx); 
[y,ye,Ny] = FDGrid(ylim(1),ylim(2),dy);

[x_eta,y_eta] = meshgrid(x,y);      % eta,H grid points
[x_u,y_u] = meshgrid(xe,y);         % u,U   grid points
[x_v,y_v] = meshgrid(x,ye);         % v,V   grid points
[x_q,y_q] = meshgrid(xe,ye);        % vorticiy grid points


betay =y_q*2.0e-11;				    % Planetary vorticity gradient (s-1)
f = f0+betay; 						% Planetary vorticity (s-1)


% other parameters
alpha = 1;                          % alpha = 1 = Include non-linear terms
 

%%% Scales of motion
L = Rd; 				    % Length scale 
U = 2.0; 					% Velocity scale (2 m/s)
T = L/U; 					% Time scale (s)

%%% Initial conditions - Anticyclonic eddies
[u0,v0,eta0] = ic(x,xe,y,ye,f0,g,300000,0,U,Rd);
u0 = u0/U;
v0 = v0/U;
eta0 = eta0*g/abs(U*L*abs(f0));

%%% Choosing obstacles geometry

if exper == 'n'
	%% No topography
	depth = 0.0;
	mask_u = ones(size(u0));
	mask_v = ones(size(v0));
	mask_q = ones(size(x_q));

	iu1 = mask_u < depth;
	iv1 = mask_v < depth;
	iq1 = mask_q < depth;

	filename = ['../matfiles/outputs_',hm,'_',exper,'topo.mat'];

else
	%% Meridional Ridge
	r = Rd/1.75;
	d = input('Distance between islands: ');
	[mask_u,mask_v,mask_q] = build_ridge(x,xe,0.0,y,ye,Ly,d,r);
	depth = 0.0;
	iu1 = mask_u == depth;
	iv1 = mask_v == depth;
	iq1 = mask_q == depth;

	filename = ['../matfiles/outputs_',hm,'_',exper,'topo_d',num2str(d/Rd,'%d'),'Rd.mat'];	
end

%%% Integration - Dimensionless form
%% Save outputs every 50 time steps
ts = [50:50:10500];

[volume,energy,enstrophy,u,v,eta,zeta] = ...
	nl_isw_2d_cb(dx/L,dy/L,xlim/L,ylim/L,f/abs(f0),alpha,H/H,g/g,dt/T,ts,u0,v0,eta0,mask_u,mask_v);

mask_u(iu1) = nan;
mask_v(iv1) = nan;
mask_q(iq1) = nan;

for t = [1:length(ts)]
	u(:,:,t) = u(:,:,t).*mask_u;
	v(:,:,t) = v(:,:,t).*mask_v;
	zeta(:,:,t) = zeta(:,:,t).*mask_q;	 
end

disp(['Saving outputs'])
eval(['save ',filename,' u0 v0 x xe y ye volume volume energy enstrophy u v eta zeta'])
