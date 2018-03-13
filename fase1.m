function [A1, b0, c0, x0, B, N] = fase1(A0,b0,c0)
% purpose: Versión del Simplex en la Fase I
% minimizar c^T x
% sujeto a Ax <= b , x >= 0 , b en R^m
%
% In : A0 ... mxn matrix
% b0 ... column vector with as many rows as A
% c0 ... column vector with as many entries as one row of A
%
% Out: A1 ... A0 modificada
% b0 ... b0 modificada
% c0 ... c0 modificada
% x0 ... una solución básica factible
% B ... Indices básicos
% N ... Indices no Básicos

%Checa cuales restricciones son no negatvas (1) o negativas (0)
if(all(b0>=0)) 
%caso trivial, pues 0 está en el Conjunto Factible

       A1=[A0,eye(size(A0,1))]; %Matriz modificada
       c1=[c0;zeros(size(A0,1),1)]; %c modificado
       b1=b0;
       x0=[zeros(size(A0,1),1); b0];
       N=[1:size(A0,2)];
       B=[(size(A0,2)+1):size(A1,2)];
       
       
else %caso complejo
    
       %restr= (b0>=0); %vector lógico
end
end

