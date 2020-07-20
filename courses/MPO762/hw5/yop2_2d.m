function g=yop2_2d(g, f, a0, a)
%
% Performs the following operation
% g(i,j-jgo) = a0*g(i,j-jgo) + a1*f(i,j) + a2*f(i,j-1) for i = 1,...,M
%                                                          j = 1,...,N
% usage is
% function g=yop2_2d(g, f, a0, a)
%                                                          j = 2,...,N+jgo
% For an y-derivative set a1=1/dy and a2=-1/dy
% For an y-average    set a1=1/2  and a2= 1/2
% Set a0=0.0_r8 to assign result to array g 
% Set a0=1.0_r8 to add result to array g 
%
% The g array offset, igo, is figured out automatically from the data
% if ubound(g,2) > ubound(f,2) then jgo=0 and we have a p-points to v-points op
%
% real(r8), intent(inout) :: g(:,:)   ! output array
% real(r8), intent(in)    :: f(:,:)   ! input array
% real(r8), intent(in)    :: a0, a(2) ! coefficients
%
%
% Author: Dr. Mohamed Iskandarani 
% Modified by Tiago Bilo (Fall 2016)


  Mg = size(g);      % size of array g
  Mf = size(f);      % size of array f

  if (Mg(2) > Mf(2)) % offset of array for the storage part
    jgo = 0;         % operation from p to v-points set jgo=0
                     % end points g untouched and must be set by BC
  else
    jgo = 1;         % operation from v to p-points set jgo=1
  end
  if (Mg(1)~=Mf(1))  % check that f and g conform in 2nd dimension
   disp('Error, 1st dimension of arrays do not conform');
  end
  N = min(Mg(2),Mf(2)); % number of x-cells
  g(:,2-jgo:N) = a0*g(:,2-jgo:N) + a(1)*f(:,2:N+jgo)+a(2)*f(:,1:N+jgo-1);
