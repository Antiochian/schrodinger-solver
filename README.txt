 CO22 - Numerical Solution of Schrodinger's Equation
---------------------------------------------------- 
 by ANTIOCHIAN 2k19 // Lets All Love Lain
---------------------------------------------------- 
Synopsis:
        CO22(delta, maximum x, order, Approximate Energy);
  
   Description:
       Inputs conditions and produces a numerical solution to the simple
       harmonic oscillator problem, which is then plotted and compared to
       the analytic solution.
   Inputs:
     -delta: the "fineness" of the numerical approximation. Smaller is
     more accurate but slower.
     -x1:    the upper bound of integration, which goes from 0 to x1
     -n:     the order of the wavefunction, integers only
     -E0:    the energy level, can be approximate

   Output:
     -none except for graph

Sample Use:  CO22 (0.05 , 5 , 1 , 3)
Sample Output:
![Sample Output](Sample_Output.png)
