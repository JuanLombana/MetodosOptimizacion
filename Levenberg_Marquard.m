clc
clear
close


syms x1 x2;
f = (x1-4)^4+(x1-4*x2)^2;
x_1(1)=0;
x_2(1)=0;
lamb=1.1;
df_dx1=diff(f,x1);
df_dx2=diff(f,x2);
d2f_dx12=diff(df_dx1,x1);
d2f_dx22=diff(df_dx2,x2);
d2f_dx1dx2=diff(df_dx1,x2);
grad = [subs(df_dx1,[x1,x2],[x_1(1),x_2(1)]), subs(df_dx2,[x1,x2],[x_1(1),x_2(1)])];
%H = [subs(d2f_dx12,[x1,x2],[x_1(1),x_2(1)]), subs(d2f_dx1dx2,[x1,x2],[x_1(1),x_2(1)]); subs(d2f_dx1dx2,[x1,x2],[x_1(1),x_2(1)]), subs(d2f_dx22,[x1,x2],[x_1(1),x_2(1)]) ];

i = 1;
while i < 3
    %l_T = inv(lamb*eye(2)+H);
    l_T = lamb^-1;
    x_1(i+1) = x_1(i)-l_T*grad';
    x_2(i+1) = x_2(i)-l_T*grad';
    z(i+1) = double(subs(f,[x1,x2],[x_1(i+1),x_2(i+1)]));
    i = i +1;
    grad = [subs(df_dx1,[x1,x2],[x_1(i),x_2(i)]), subs(df_dx2,[x1,x2],[x_1(i),x_2(i)])];
    %H = [subs(d2f_dx12,[x1,x2],[x_1(i),x_2(i)]), subs(d2f_dx1dx2,[x1,x2],[x_1(i),x_2(i)]); subs(d2f_dx1dx2,[x1,x2],[x_1(i),x_2(i)]), subs(d2f_dx22,[x1,x2],[x_1(i),x_2(i)]) ];
end

fcontour(f, 'Fill', 'On');
hold on;
plot(x_1,x_2,'*-r');
grid on;

X1 = x_1';
X2 = x_2';
Iterations = (1:i)';
Z = z';
T = table(Iterations,X1,X2,Z);
fprintf('Minimo: %f\n\n',subs(f,[x1,x2],[x_1(i),x_2(i)]));
disp(T)