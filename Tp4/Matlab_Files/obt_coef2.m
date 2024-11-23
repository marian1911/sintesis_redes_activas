function [k,W0,Q] = obt_coef2(a,denom)

%obt_coef
%Esta función toma el numerador y el denominador de una función de transferencia de segundo orden y devuelve la ganancia "k", el "W0" y el “Q” del correspondiente sistema prototipo de segundo orden
%
%usted debe llamar a la funcion obt_coef(a,denom) y poner en lugar de 'a'
%el numerador de la funcion de segundo orden y en el lugar de 'denom' el
%vector que representa el denominador de dicha funcion de segundo orden
% 
% see also STEP

num = a/denom(1);
denom = denom/denom(1);
W0 = sqrt(denom(3));
k = num(1);
z = denom(2)/(W0);
Q = 1/z;

end