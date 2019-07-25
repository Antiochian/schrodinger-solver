function [ finalvalue ] = solve_eigenvalue (E0,n)
% EIGENVALUE SOLVER
%---------------------------------------------------- 
% by HEIRON 2k19 // Lets All Love Lain
%----------------------------------------------------
%   Synopsis:
%        solve_eigenvalue(Approx Energy, n);
%  
%   Description:
%       Inputs a guess for an energy eigenvalue and outputs the closest
%       eigenvalue for the odd and even case
%   Inputs:
%     -E0: best guess for E (energy eigenvalue)
%     -n : n-value for desired eigenvalue (optional)
%
%   Output:
%     -finalvalue: Revised value for E
%

%-----------------------
%SET UP OF KEY VARIABLES
%-----------------------
%[ sample variables suggested in script ]
xarray=linspace(0,5,101);

x0= xarray(1);     %Extract the first and last values from the given xarray
xf = xarray(end);
delta = (xf-x0)/size(xarray,2); %Extract the spacing
jmax = size(xarray,2)-1;        %Set up index (note the -1 index correction)
step = delta;                   %intial step (in +ve direction)
accuracy = 0.01;                %choose accuracy desired
%choose direction of step (toward closest integer)
if mod(E0,1) < 0.5 && E0>0.5
    step = -1*step;
else
end

%n detection
if nargin == 1 %detect if n is inputted or not
    psi0 = 0; %DEFAULT to odd n
    dpsi0 = 1;
else
    if mod(n,2) == 0; %i.e: if n is even
    psi0 = 1;
    dpsi0 = 0;
    else
    psi0 = 0;
    dpsi0 = 1;
    end
end

E = E0; %Set starting value of E to E0

%------------------------
%CALCULUS + TAYLOR APPROX.
%------------------------
%This is done with zero prior knowledge in order to ensure program is as
%general as possible and doesnt just work for Schrodinger Eqs.
%so we find df/dx terms using symbolic variables:
f = @(x) x^2-E;
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

%NUMEROV CALCULATION
%treat as A*Psi(+1)=B*Psi(0)+C*Psi(-1) to simplify + aid in debugging
PrevFinalPsi = 0; %we choose 0 as a starting value to ensure FV/Psi !< 0
go = true; 
while go
PsiVector = solve_numerov ( f , xarray, psi0, dpsi0 );
FinalPsi = PsiVector(end); %take only value of wavefunction at x=xf

if PrevFinalPsi/FinalPsi < 0 %CHECK FOR SIGN CHANGE!
        if abs(FinalPsi) < accuracy;
            break %stop when as close as required
        else
            step = -0.5*step; %increase finetuning granularity, REVERSE DIRECTION
        end
elseif PrevFinalPsi/FinalPsi < 1.2 && abs(PrevFinalPsi/FinalPsi) > 0 %if value is barely changing!
    step = 2*step;
else %else carry on as usual
end



%update values for next step
E = E + step;
PrevFinalPsi = FinalPsi; %set new FinalValue
f = @(x) x^2-E;

end
%--------
%OUTPUT
%--------
finalvalue = E; %output Eigenvalue

end