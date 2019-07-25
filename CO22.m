function [ output ] = CO22 (delta , x1 , n , E0 )
% CO22 - Numerical Solution of Schrodinger's Equation
%---------------------------------------------------- 
% by HEIRON 2k19 // Lets All Love Lain
%---------------------------------------------------- 
%Synopsis:
%        CO22(delta, maximum x, order, Approximate Energy);
%  
%   Description:
%       Inputs conditions and produces a numerical solution to the simple
%       harmonic oscillator problem, which is then plotted and compared to
%       the analytic solution.
%   Inputs:
%     -delta: the "fineness" of the numerical approximation. Smaller is
%     more accurate but slower.
%     -x1:    the upper bound of integration, which goes from 0 to x1
%     -n:     the order of the wavefunction, integers only
%     -E0:    the energy level, can be approximate
%
%   Output:
%     -none except graph
%
%   Sample Use:  CO22 (0.05 , 5 , 1 , 3)

%-------------------
%CORRECT EIGENVALUE:
%-------------------
E=solve_eigenvalue (E0,n);
%-----------------------
%SET UP FOR FUNCTIONCALL
%-----------------------
f = @(x) x^2-E; %set up function (this is for Schrodinger only)
jmax = x1/delta; %count number of iterations required
x = linspace(0, x1, jmax+1); %create x-vector

%Set initial conditions
if mod(n,2) == 0; %i.e: if n is even
    psi0 = 1;
    dpsi0 = 0;
else
    psi0 = 0;
    dpsi0 = 1;
end

%---------------------
%NUMEROV FUNCTIONCALL
%---------------------

psi = solve_numerov(f, x, psi0, dpsi0); %Functioncall!

%-----------------
%ANALYTIC SOLUTION
%-----------------

AnalyticSol=zeros(jmax+1,1); %Avoid concatenation to increase speed
    for i=0:jmax
        AnalyticSol(i+1,1) = hermiteH(n,i*delta)*exp(-((i*delta).^2)/2);
    end
    
%note that there will be a corrective factor required
%to find this we take the average ratio of points and find the MODAL average:
error = psi./AnalyticSol; %note this is a vector of errors
errorratio = round(error,2,'significant'); %Round to 2 sig. figs (changeable)
scale = mode(errorratio); %We take the mode because the error hugely spikes near xf and warps the mean
%-----------
%PLOTTING
%-----------
plot(x, psi); %plot numerical solution (blue)
hold on;
plot(x,scale*AnalyticSol(:,1),':r'); %plot analytic solution (red, dashed)
hold off;
xlabel('x'),ylabel('Psi(x)'),axis square, grid on
end