% Author: Tiago Bilo
% CFD - Fall 2016
% Class project:
% Self-propagating eddies and obstacles interactions in a zonal channel 
% solved by the Inviscid Shallow Water equations 
%
% the idea is to observe what kind of behavior we can recover from a simple inviscid model
% when compared to lab experiments and NBC rings observations (Cenedese et al., 2005) !

% Clear memory and close all figures
clear all 
close all 



%%%% Loading topographic mask
load '../matfiles/topo.mat'
load '../matfiles/mask_model_cood.mat'


disp(['Setting model parameters'])
%%%% Model parameters 
% Physical parameters
H = 1000.0;                       % Cte depth (1000 m)
g = 0.0045;                       % Reduced gravity gravity (<< 10 m/s2 and < 0) or acceleration due to gravity (~10 m/s2)
f0 = gsw_f(13);		  		  % Coriolis parameter (s-1)
Rd = sqrt(H*g)/f0;				  % Rossby Radius of deformation

Lx = 1000000.0;				      % Channel length (m) 
Ly = 1200000.0;                    % Channel width (m)


% I have to figure out the dimensions of the islands
% and the distance between them 

% Grid parameters (Arakawa C-grid)
dx = 5000;            		  			        % zonal resolution 5 km (m)
dy = dx;                          				% meridional resolution (m)
dt = 900;                       	            % Time step of 15 min (s) 

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

%%% Initial conditions - Dimensionless form
[u0,v0,eta0] = ic(x,xe,y,ye,f0,g,100000,-250000,U,Rd);
u0 = u0/U;
v0 = v0/U;
eta0 = eta0*g/abs(U*L*f0);


%%% Creating the topographic masks
mask_u = griddata(double(x_mask), double(y_mask),double(topo), x_u, y_u); 
mask_v = griddata(double(x_mask), double(y_mask),double(topo), x_v, y_v);

depth = 1000;

iu1 = mask_u < depth;
iu2 = mask_u >= depth;

mask_u(iu1) = 0.0;
mask_u(iu2) = 1.0;

iv1 = mask_v < depth;
iv2 = mask_v >= depth;

mask_v(iv1) = 0.0;
mask_v(iv2) = 1.0;

disp(['Integrating model'])
%%% Integration - Dimensionless form
%% Save outputs every 100 time steps
ts = [100:100:6000];

[volume,energy,enstrophy,u,v,eta] = ...
	nl_isw_2d_cb(dx/L,dy/L,xlim/L,ylim/L,f/abs(f0),alpha,H/H,g/g,dt/T,ts,u0,v0,eta0,mask_u,mask_v);


mask_u(iu1) = nan;
mask_v(iv1) = nan;

for t = [1:length(ts)]
	u(:,:,t) = u(:,:,t).*mask_u;
	v(:,:,t) = v(:,:,t).*mask_v; 
end
