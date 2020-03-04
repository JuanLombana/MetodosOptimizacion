clc
clear
close

syms x1 x2 x3;
fo = x1^2+8*x2^2+3*x1*x3;
f = (2*x2-x3-1)^2+8*x2^2+3*(2*x2-x3-1)*x3;
x =[0,0;];
df_dx1=diff(f,x2);
df_dx2=diff(f,x3);
d2f_dx12=diff(df_dx1,x2);
d2f_dx22=diff(df_dx2,x3);
d2f_dx1dx2=diff(df_dx1,x3);
grad = [subs(df_dx1,[x2,x3],x(1,:)), subs(df_dx2,[x2,x3],x(1,:))];
H = [subs(d2f_dx12,[x2,x3],x(1,:)), subs(d2f_dx1dx2,[x2,x3],x(1,:)); subs(d2f_dx1dx2,[x2,x3],x(1,:)), subs(d2f_dx22,[x2,x3],x(1,:)) ];

i = 1;
while i < 3
    x(i+1,:) =(x(i,:)'-H\grad')';
    z(i+1) = double(subs(f,[x2,x3],x(i+1,:)));
    i = i +1;
    grad = [subs(df_dx1,[x2,x3],x(i,:)), subs(df_dx2,[x2,x3],x(i,:))];
    H = [subs(d2f_dx12,[x2,x3],x(i,:)), subs(d2f_dx1dx2,[x2,x3],x(i,:)); subs(d2f_dx1dx2,[x2,x3],x(i,:)), subs(d2f_dx22,[x2,x3],x(i,:)) ];
end

fcontour(f, 'Fill', 'On');
hold on;
plot(x(:,1), x(:,2),'*-r');
grid on;
X1 = 2*x(:,1)-x(:,2)-1;
X2 = x(:,1);
X3 = x(:,2);
Iterations = (1:i)';
Z = z';
T = table(Iterations,X1,X2,X3,Z);
fprintf('Minimo: %f\n\n',subs(f,[x2,x3],x(i,:)));
fprintf('Minimo Original: %f\n\n',subs(fo,[x1,x2,x3],[X1(i),X2(i),X3(i)]));
disp(T)

