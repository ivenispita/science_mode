function p = rk3(p,c,dx,dt,periodic)
%
% time steps the variable h through a single step of size dt
% usage is
% p = rk3(p,c,dx,dt,periodic)
% where p: is the variable at the present time level
%          and is the variable at the next time level
%       dt: is the time step
%       c: is the constant advection speed
%       dx: is the grid size
%       periodic: is the flag for periodic BC
%

pt = zeros(size(p));       % temporary mid-stage solution

% Stage 1
r = rhsadv(p,c,dx,periodic); % tendency for stage 1
pt = p + dt*r;               % stage 1 update

% Stage 2
r = rhsadv(pt,c,dx,periodic);   % tendency for stage 2
pt = 0.75*p + 0.25*(pt + dt*r); % Stage 2 update

% Stage 3
r = rhsadv(pt,c,dx,periodic);   % tendency for stage 3
p = (p + 2.0*(pt + dt*r))/3.0;  % final update
