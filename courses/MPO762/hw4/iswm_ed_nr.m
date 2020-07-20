% Clear memory and close all figures
clear all 
close all 


% Model parameters
H1 = 10.0;                        % Cte depth [m]
a = 1280.0;                       % Basin x-dimension [m]
b = 640.0;                        % Basin y-dimension [m]
alpha = 0.1;                      % Wave amplitude [m]
l2 = 25.0*25.0;                   % Wave height decay scale [m2]

g = 10;                           % Acceleration due to gravity [m/s2]

dt = 0.1;                         % Time step [s]                     
dx = [2];                     % Zonal grid spacing [m]
dy = dx;                          % Meridional grid spacing [m]

ts = [32,64,128,256];             % Save solutions time-steps [s] 

% Run the 2D-Inviscid Shallow Water Model
% in a closed rectangular basin described by a C-grid
% 
% The model integration is done by calling isw_2d_cp.m
for i = [1:length(dx)]

    % Grid cells and grid cell edges - Grid set-up
    [x,xe,Nx] = FDGrid(0,a,dx(i)); 
    [y,ye,Ny] = FDGrid(0,b,dy(i));

    [x_eta,y_eta] = meshgrid(x,y);      % eta,H grid points
    [x_u,y_u] = meshgrid(xe,y);         % u,U   grid points
    [x_v,y_v] = meshgrid(x,ye);         % v,V   grid points
    [xe,ye] = meshgrid(xe,ye);          % vorticiy grid points
    
    
    % Variable topography
    H2 = 10.0+5.0*tanh((x_eta-100.0)/16.0);
    
    % Square of the radial distance from the center of the basin
    r2 = (((a/2)-x_eta).^2)+(((b/2)-y_eta).^2);
     
    % Printing some Grid information
    disp(' ')    
    disp(' ')
    disp(['[outputs ',num2str(i),']'])
    disp(['Resolution:      ',num2str(dx(i)),' m'])
    disp(['Number of x-cells: ',num2str(Nx)])
    disp(['Number of y-cells: ',num2str(Ny)])    
    disp(' ')


    % Initial condition
    u0 = zeros(size(x_u));
    v0 = zeros(size(x_v));
    eta0 = alpha.*exp(-r2/l2);


    %% Integration
    % Constant depth
    [volume1,energy1,u1,v1,eta1] = ...
        isw_2d_cb(dx(i),dy(i),a,b,H1,g,dt,ts,u0,v0,eta0);


    [volume2,energy2,u2,v2,eta2] = ...
        isw_2d_cb(dx(i),dy(i),a,b,H2,g,dt,ts,u0,v0,eta0);


    % Reporting results
    disp('Constant Depth H1 = 10 m')
    disp('_________________________________________________')
    disp(['Time step (s):           ',num2str(ts)])
    disp('_________________________________________________')
    disp(['Volume displacement (m^3/s):   ',num2str(volume1)])    
    disp(['Total energy (kJ/density unit): ',num2str(energy1/1000)])
    disp(' ')
    
    disp('Variable Depth H2')
    disp('_________________________________________________')
    disp(['Time step (s):           ',num2str(ts)])
    disp('_________________________________________________')
    disp(['Volume displacement (m^3/s):   ',num2str(volume2)])    
    disp(['Total energy (kJ/density unit): ',num2str(energy2/1000)])
    disp(' ')   


end