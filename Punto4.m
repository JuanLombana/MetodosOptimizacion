%PRUEBE QUE EL GRADIENTE DE c^Tx = c

%definir las variables simbolicas
syms x1 x2 x3 x4;
%función objetivo
f = -1*x1+0*x2+7*x3+4*x4;
%calcular el gradiente de la función
g = ([ diff(f,x1), diff(f,x2), diff(f,x3),diff(f,x4)])';
%valor esperado
c = [-1;0;7;4];
%mostrar la comparación
disp(g==c);

%PRUEBE QUE EL JACOBIANO DE (Ax-b) = A^T
%definir la función restricción (Ax-b)
h1 = 1*x1+1*x2+1*x3+1*x4-4;
%calcular el jacobiano de la restricción
j_h = [diff(h1,x1); diff(h1,x2); diff(h1,x3); diff(h1,x4)];
%Valor esperado
AT=([1,1,1,1])';
%mostra la compación
disp(AT== j_h);
