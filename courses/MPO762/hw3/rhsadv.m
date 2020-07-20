function r = rhsadv(p,c,dx,periodic)
%
% Compute time-tendency of finite difference advection term of the form, p_t = -c p_x
% usage is
% r = rhsadv(p,c,dx,periodic)
% where p is the variable at the present time level
%       u is the flow velocity at cell edges
%       dx is the grid size
%       periodic is a flag for periodic BC.
%
% Created by Mohamed Iskandarani
% Modified by Tiago Bilo (Fall 2016)


% r = cd2(p,dx, periodic);          % 1st derivative, centered 2nd order
% r = cd4(p,dx, periodic);          % 1st derivative, centered 4th order
% r = cd6(p,dx, periodic);          % 1st derivative, centered 6th order
% r = cd8(p,dx, periodic);          % 1st derivative, centered 8th order
  r = bd5(p,dx, periodic);          % 1st derivative, backward 5th order


  r =-c*r;                          % advection term
