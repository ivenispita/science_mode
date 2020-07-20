function du = bd5(ue,dx,periodic)
%
% Returns the 5th-order backward difference derivative of the function 
% defined in the array ue
% Usage is
% du = bd5(ue,dx,periodic)
% where 
%  dx        is the grid spacing (assumed uniform)
%  periodic  is a flag which is set to 1 to indicate periodic treatment of ghost values
%  ue        is the array to be differentiated 
%  du        is the array containing the finite difference approximation
%
  Ne=length(ue);
  du(4:Ne-2) =(-2*ue(1:Ne-5)+15*ue(2:Ne-4)...
               -60*ue(3:Ne-3)+20*ue(4:Ne-2)...
               +30*ue(5:Ne-1)-3*ue(6:Ne))/(60*dx); % backward 5th order


  if (periodic)
    du(3) = (-2*ue(Ne-1)+15*ue(1)...
             -60*ue(2)+20*ue(3)...
             +30*ue(4)-3*ue(5))/(60*dx); % backward 5th order

    du(2) = (-2*ue(Ne-2)+15*ue(Ne-1)...
             -60*ue(1)+20*ue(2)...
             +30*ue(3)-3*ue(4))/(60*dx); % backward 5th order

    du(1) = (-2*ue(Ne-3)+15*ue(Ne-2)...
             -60*ue(Ne-1)+20*ue(1)...
             +30*ue(2)-3*ue(3))/(60*dx); % backward 5th order


    du(Ne-1) = (-2*ue(Ne-4)+15*ue(Ne-3)...
             -60*ue(Ne-2)+20*ue(Ne-1)...
             +30*ue(Ne)-3*ue(2))/(60*dx); % backward 5th order

    du(Ne) = (-2*ue(Ne-3)+15*ue(Ne-2)...
             -60*ue(Ne-1)+20*ue(Ne)...
             +30*ue(2)-3*ue(3))/(60*dx); % backward 5th order

  end 

