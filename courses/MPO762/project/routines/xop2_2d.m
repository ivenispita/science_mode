function g=xop2_2d(g, f, a0, a)
%
% Performs the following operation on a C-grid
% g(i-igo,j) = a0*g(i-igo,j) + a1*f(i,j) + a2*f(i-1,j) for i = 2,...,M+igo
%                                                          j = 1,...,N
% usage is
% function g=xop2_2d(g, f, a0, a)
%
% For an x-derivative set a1=1/dx and a2=-1/dx
% For an x-average    set a1=1/2  and a2= 1/2
% Set a0=0._r8 to assign result to array g 
% Set a0=1._r8 to add result to array g 
% The g array offset, igo, is figured out automatically from the data
% if ubound(g,1) > ubound(f,1) then igo=0 p-points to u-points op
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

  if (Mg(1) > Mf(1)) % offset of array g
    igo = 0;         % operation from p to u-points, set igo=0
                     % end points g untouched and must be set by BC
  else
    igo = 1;         % operation from u to p-points set, igo=1
  end
  if (Mg(2)~=Mf(2))  % check that f and g conform in 2nd dimension
   disp('Error, 2nd dimension of arrays do not conform');
  end
  M = min(Mg(1),Mf(1)); % number of x-cells
  g(2-igo:M,1:Mg(2)) = a0*g(2-igo:M,1:Mg(2)) + a(1)*f(2:M+igo,1:Mg(2)) + a(2)*f(1:M+igo-1,1:Mg(2));
