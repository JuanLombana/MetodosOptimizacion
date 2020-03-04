clc
clear
close

syms x1 x2;
f = (x1-4)^4+(x1-4*x2)^2;
x=[0,0;1/10000000000000,0;];
df_dx1=diff(f,x1);
df_dx2=diff(f,x2);
grad_a = [subs(df_dx1,[x1,x2],x(1,:)), subs(df_dx2,[x1,x2],x(1,:))];
grad = [subs(df_dx1,[x1,x2],x(2,:)), subs(df_dx2,[x1,x2],x(2,:))];
i = 2;
D=eye(2);
while i < 4
    s = (x(i,:)-x(i-1,:))';
    q = grad'-grad_a';
    d_new = D + ((s*s')/(q'*s))-((D*(q*q')*D')/(q'*D*q));
    x(i+1,:) = x(i)-d_new*grad';
    z(i+1) = double(subs(f,[x1,x2],x(i+1,:)));
    i = i +1;
    grad_a = grad;
    grad = [subs(df_dx1,[x1,x2],x(i,:)), subs(df_dx2,[x1,x2],x(i,:))];
    D = d_new;
end

fcontour(f, 'Fill', 'On');
hold on;
plot(x(:,1), x(:,2),'*-r');
grid on;

X1 = x(:,1);
X2 = x(:,2);
Iterations = (1:i)';
Z = z';
T = table(Iterations,X1,X2,Z);
fprintf('Minimo: %f\n\n',subs(f,[x1,x2],x(i,:)));
disp(T)