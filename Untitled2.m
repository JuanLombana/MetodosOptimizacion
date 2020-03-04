A=[1 1];
b=6;
x0=[0 0]';
options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); 
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(@funZ,x0,A,b,[],[],[],[],@funZ_nonlcon,options);
disp("mínimo: " + fval)
disp("X:");S
disp(x);
