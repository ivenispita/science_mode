function [x,xe,dx] = FDGrid(xmin,xmax,N)
    %
    % Returns the coordinates of the cell (x), grid cell edges (xe), and the spacing (dx) 
    % of an equally spaced 1-D numerical grid.  
    % 
    % 
    % Input:  
    %    xmin and xmax: interval limits of the grid 
    %         
    %    N: Number of grid cells 
    %   
    % Example: 
    %
    %            >> [x,xe,dx] = FDGrid(0,10,5)
    %            x = 
    %                [1,3,5,7,9]
    %             
    %            xe = 
    %                [0,2,4,6,8,10]
    %
    %      
    %            dx =             
    %                 2
    %
    %    
    %    Grid: 
    %    
    %   0       1     2     3      4      5     6     7     8      9    10
    %    |       .     |     .      |      .     |     .     |      .     |
    %  xe1     x1    xe2    x2    xe3    x3    xe4   x4    xe5     x5    xe6         
    %
    % 
    % Author: Tiago Bilo
    % CFD - Fall 2016 
    % Problem set 1: 
    % 2. Programing discretization of a grid 

    % Grid spacing
    dx = (xmax-xmin)/N;

    % Grid cell coordinates
    x = [xmin+(dx/2.0):dx:xmax-(dx/2.0)];

    % Coordinates of the grid cell edges 
    xe = [xmin:dx:xmax];



