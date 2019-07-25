function [ psi ] = solve_numerov (f , xarray , psi0 , dpsi0 )
% NUMEROV PLOTTER by HEIRON 2k19 // Lets All Love Lain
%   Synopsis:
%        solve_numerov(f(x), x, psi(0), dpsi(0));
%  
%   Description:
%       Inputs second-order linear ODEs with no first-order terms and
%       outputs a numerical solution as a vector.
%   Inputs:
%     -f: The function to be called on the right hand side of the equation (see hint
%      section), receives x and returns the value of f(x).
%     -xarray: array of the integration, the first element is the lower bound, x0, and the last
%      element is the upper bound, x1.
%     -psi0: the value of \psi at x0, i.e. \psi(x0)
%     -dpsi0: the value of \psi derivative at x0, i.e. d\psi(x0)/dx
%
%   Output:
%     -psi: array of values of \psi with each element in the array corresponding to the same
%      element in x
%


%We are integrating from x0 to xf with a spacing of delta
%-----------------------
%SET UP OF KEY VARIABLES
%-----------------------
x0= xarray(1);     %Extract the first and last values from the given xarray
xf = xarray(end);
delta = (xf-x0)/size(xarray,2); %Extract the spacing
jmax = size(xarray,2)-1;        %Set up index (note the -1 index correction)


%------------------------
%CALCULUS + TAYLOR APPROX.
%------------------------
%This is done with zero prior knowledge in order to ensure program is as
%general as possible and doesnt just work for the Schrodinger Eq.
%so we manually find df/dx terms using symbolic variables:
syms x;
syms y;
syms z;
d1= diff(f(x));       %symbolic equation for f'(x)
d1f= @(y) subs(d1,y); %define f'(x) as lambda function
d2 = diff(d1f(y));    %symbolic equation for f''(x)
d2f= @(z) subs(d2,z); %define f''(x) as lambda function

%Use taylor to find initial Psi conditions
%NB: This is split up only to make debugging/transcribing easier
TaylorEven = psi0+(delta.^2)*f(0)*psi0/2+(delta.^4)*(d2f(0)*psi0+2*d1f(0)*dpsi0+((f(0)).^2)*psi0 )/24;
TaylorOdd  = (delta*dpsi0 + (delta.^3)*(f(0)*dpsi0+d1f(0)*psi0) )/6;

%Sum the two parts together for our approximation of psi(delta)
psi1 = double(TaylorEven + TaylorOdd);

%-----------------------
%ITERATIVE STEP
%-----------------------
%create 2-wide matrix NumericalSol to store Psi(j) and f(j)
NumericalSol = zeros(jmax+1,2); %We create a zero matrix instead of concatenating for speed reasons

%Input initial conditions for Psi [stored in col 1]
NumericalSol(1,1) = psi0;
NumericalSol (2,1) = psi1;
%Input f(x) for each value of x given into same matrix
for i=0:jmax
NumericalSol(i+1,2) = f(i*delta);
end

%NUMEROV CALCULATION
%treat as A*Psi(+1)=B*Psi(0)+C*Psi(-1) to simplify + aid in debugging
for j=2:jmax
    A =  (1 - (delta.^2)*NumericalSol(j+1,2)/12);
    B =  (2 + (delta.^2)*NumericalSol(j ,2)*5/6);
    C = -(1 - (delta.^2)*NumericalSol(j-1,2)/12);
    %Actual Numerov iteration step:
   NumericalSol(j+1,1)= (B*NumericalSol(j,1) + C*NumericalSol(j-1,1))/A;
end
%-----------------------
%OUTPUT
%-----------------------
psi = NumericalSol(:,1); %output array as desired

end