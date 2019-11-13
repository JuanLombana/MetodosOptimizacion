clc
clear
close

syms x1 x2;
f = (x1-4)^4+(x1-4*x2)^2;
x_1(1)=0;
x_1(2)=0.001;
x_2(1)=0;
x_2(2)=0.001;
df_dx1=diff(f,x1);
df_dx2=diff(f,x2);
grad_a = [ 0, 0];
grad = [subs(df_dx1,[x1,x2],[x_1(1),x_2(1)]), subs(df_dx2,[x1,x2],[x_1(1),x_2(1)])];
i = 2;
H=eye(2);
while i < 4   
    s = [x_1(i)-x_1(i-1);x_2(i)-x_2(i-1)];
    q = grad'-grad_a';
    h_new = H + ((q*q')/(q'*s))-((H*(s*s')*H')/(s'*H*s));
    inv_h = inv(h_new);
    x_1(i+1) = x_1(i)-inv_h(1,:)*grad';
    x_2(i+1) = x_2(i)-inv_h(2,:)*grad';
    z(i+1) = double(subs(f,[x1,x2],[x_1(i+1),x_2(i+1)]));
    i = i +1;
    grad_a = grad;
    grad = [subs(df_dx1,[x1,x2],[x_1(i),x_2(i)]), subs(df_dx2,[x1,x2],[x_1(i),x_2(i)])];
    H = h_new;
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