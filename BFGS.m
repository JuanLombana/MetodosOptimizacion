clc
clear
close

syms x1 x2;
f = (x1-4)^4+(x1-4*x2)^2;
x_1(1)=0;
x_2(1)=0;
df_dx1=diff(f,x1);
df_dx2=diff(f,x2);
grad_a = [ 0, 0];
grad = [subs(df_dx1,[x1,x2],[x_1(1),x_2(1)]), subs(df_dx2,[x1,x2],[x_1(1),x_2(1)])];
i = 1;
H_1=?;
H_2=?;
while i < 3
    s_1 = x_1(i)-x_1(i-1);
    s_2 = x_2(i)-x_2(i-1);
    q = grad'-grad_a';
    h_new_1 = H_1 + (q*q'/q'*s_1)-(H_1*(s_1*s_1')*H_1'/s_1'*H_1*s_1);
    h_new_2 = H_2 + (q*q'/q'*s_2)-(H_2*(s_2*s_2')*H_2'/s_2'*H_2*s_2);
    x_1(i+1) = x_1(i)-inv(h_new_1)*grad';
    x_2(i+1) = x_2(i)-inv(h_new_2)*grad';
    z(i+1) = double(subs(f,[x1,x2],[x_1(i+1),x_2(i+1)]));
    i = i +1;
    grad_a = grad;
    grad = [subs(df_dx1,[x1,x2],[x_1(i),x_2(i)]), subs(df_dx2,[x1,x2],[x_1(i),x_2(i)])];
    H_1 = h_new_1;
    H_2 = h_new_2;
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