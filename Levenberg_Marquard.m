clc
clear
close

syms x1 x2;
f = (x1-4)^4+(x1-4*x2)^2;
x=[0,0;];
lamb=1.1;
df_dx1=diff(f,x1);
df_dx2=diff(f,x2);
d2f_dx12=diff(df_dx1,x1);
d2f_dx22=diff(df_dx2,x2);
d2f_dx1dx2=diff(df_dx1,x2);
grad = [subs(df_dx1,[x1,x2],x(1,:)), subs(df_dx2,[x1,x2],x(1,:))];
H = [subs(d2f_dx12,[x1,x2],x(1,:)), subs(d2f_dx1dx2,[x1,x2],x(1,:)); subs(d2f_dx1dx2,[x1,x2],x(1,:)), subs(d2f_dx22,[x1,x2],x(1,:)) ];

i = 1;
while i < 3
    h_new = inv((lamb*eye(2))+H);
    x(i+1,:) = (x(i,:)'-h_new*grad')';
    z(i+1) = double(subs(f,[x1,x2], x(i+1,:)));
    i = i +1;
    grad = [subs(df_dx1,[x1,x2], x(i,:)), subs(df_dx2,[x1,x2], x(i,:))];
    H = [subs(d2f_dx12,[x1,x2], x(i,:)), subs(d2f_dx1dx2,[x1,x2], x(i,:)); subs(d2f_dx1dx2,[x1,x2], x(i,:)), subs(d2f_dx22,[x1,x2], x(i,:)) ];
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