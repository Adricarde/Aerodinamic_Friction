%Resolucion de la EDO para las condiciones de contorno, alla vamos
clc
clear all
close all
load Vh
Uh(1,T/incr_t)=0;Uh(1,1)=H;
Uh(1,2)=Uh(1,1)-incr_t*Ve; V(1,1)=-fvel(Uh,1,Ve,H); V(1,2)=-fvel(Uh,2,Ve,H);
[V(1,1),q]=fvel(Uh,1,Ve,H);V(1,1)=-V(1,1);qt(:,1)=q;
[V(1,2),q]=fvel(Uh,2,Ve,H);V(1,2)=-V(1,2);qt(:,2)=q;
for i=2:T/incr_t;
    Uh(1,i+1)=Uh(1,i)+incr_t*0.5*(3*fvel(Uh,i,Ve,H)-fvel(Uh,i-1,Ve,H));
    Uh(1,i+1)=Uh(1,i)+incr_t*0.5*(fvel(Uh,i,Ve,H)+fvel(Uh,i+1,Ve,H));
    [V(1,i+1),q]=fvel(Uh,i+1,Ve,H);V(1,i+1)=-V(1,i+1);qt(:,i+1)=q;
end
t=0:incr_t:T;
%figure(1)
%plot(t,Uh),grid on,axis tight;
%figure(2)
%plot(t,V),grid on;
%load densidades
%qt=[qt;flipud(qt)];
%plot(qt),grid on, axis tight;
clear q, clear i;L=T/incr_t;
save Calor
disp('Condiciones de contorno calculadas')
run Principal