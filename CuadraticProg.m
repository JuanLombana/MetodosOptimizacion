H = [ 2,0; 0 8];
c = [-8, -16];
A = [1,1;1,0;-1,0;0,-1];
b = [5,3,0,0];
[x,fval]=quadprog(H,c,A,b);
