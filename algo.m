c=[-1;-1];
A=[1,2;1, 1/2;-1,-1];
b=[10;5;-1];
[x0, z0, ban, iter]= mSimplexRobusto(A, b, c)
