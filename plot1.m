clc
clear
tu=load('1.txt');
A=tu(:,1);
B=tu(:,2);
C=tu(:,3);
D=tu(:,4);
figure(4)
semilogy(A,B,'b:o');
semilogy(C,D,'b:o');