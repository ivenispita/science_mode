function du = cd4(ue,dx,periodic)
%
% Returns the 4th-order centered difference derivative of the function 
% defined in the array ue
% Usage is
% du = cd4(ue,dx,periodic)
% where 
%  dx        is the grid spacing (assumed uniform)
%  ue        is the array to be differentiated 
%  periodic  is a flag which is set to 1 to indicate periodic treatment of ghost values
%  du        is the array containing the finite difference approximation
%
  Ne =length(ue);  % number of edges
  du(3:Ne-2) = ( 8*(ue(4:Ne-1)-ue(2:Ne-3)) - (ue(5:Ne  )-ue(1:Ne-4)) ) /(12*dx); % Centered 4th order
  if (periodic==1)
    du(2   ) = ( 8*(ue(3   )-ue(1   )) - (ue(4   )-ue(Ne-1)) ) /(12*dx); % Centered 4th order
    du(1   ) = ( 8*(ue(2   )-ue(Ne-1)) - (ue(3   )-ue(Ne-2)) ) /(12*dx); % Centered 4th order
    du(Ne-1) = ( 8*(ue(Ne  )-ue(Ne-2)) - (ue(2   )-ue(Ne-3)) ) /(12*dx); % Centered 4th order
    du(Ne  ) = ( 8*(ue(2   )-ue(Ne-1)) - (ue(3   )-ue(Ne-2)) ) /(12*dx); % Centered 4th order
  end
