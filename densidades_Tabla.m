clc
clear all
close all
%f1=inline('sin(pi*x)','x');
f2=inline('1/(1+25*x*x)','x');cn=0;inc_barrido=1;


load Vh
N=size(Vh,2);
q=3;


[If2,x2,fr2,Rn2,xj2]=lagrange_equiespaciados_gradoq(N,inc_barrido,q,Vh);


figure(1)
loglog(nvect,Rmax(1,:),nvect,Rmax(2,:),nvect,Rmax(3,:),nvect,Rmax(4,:)),grid on, axis tight;
xlabel('Ln(N)');ylabel('Ln(Rn)max');legend('q=1','q=2','q=3','q=N');
figure(2)
loglog(nvect,Pimax(1,:),nvect,Pimax(2,:),nvect,Pimax(3,:),nvect,Pimax(4,:)),grid on, axis tight;
xlabel('Ln(N)');ylabel('Ln(Pi)max');legend('q=1','q=2','q=3','q=N');