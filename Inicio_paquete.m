clc
clear all
close all
%Datos básicos
H=125000;Ve=7500;R=0.15;Sw=2*pi*R*R;%Metros, Metros/Segundo, Metros, Metros*Metros, 
lambda_entrada=-pi/2;%Radianes
H=input('Altura de reentrada');Ve=input('Velocidad inicial de reentrada');R=input('Radio de la esfera');
lambda_entrada=input('Angulo de reentrada(negativo siempre)');
%Anular las dos líneas anteriores si no se desean modificar los datos
%preestablecidos


%Intervalo total de tiempo y su incremento
T=30;incr_t=0.05;
%Discretización del radio, valor del radio, discretización del ángulo y
%valor del intervalo de fi
N=10;a=R;M=10;b=pi*(1+0.55/M);
N=input('Valor de N(radial)');M=input('Valor de M(angular)');%Anular esta linea si no se desea modificar estos datos
ro=2700;kterm=230;cal_esp=900;%Densidad, conductividad termica y calor específico
c=(kterm/(ro*cal_esp))^0.5;
Tamb=300;%Temperatura inicial
mu=1.9*10^-5;%viscosidad dinámica


atmosfera_isa=[11000*ones(1,T/incr_t+1) ; 20000*ones(1,T/incr_t+1);32000*ones(1,T/incr_t+1);47000*ones(1,T/incr_t+1);51000*ones(1,T/incr_t+1);71000*ones(1,T/incr_t+1);84000*ones(1,T/incr_t+1)];
%Empezamos calculando la velocidad en función de la altura
load densidades %Estas densidades se pueden modificar sin ningun problema si se presentan igual que las de defecto
dens=dens(1:H/100+1,1:2);
h=dens(1:H/100+1,1);h=h';dens(:,2)=dens(:,2)*1000;
for i=1:H/100;
    drodh=(dens(i+1,2)-dens(i,2))/(dens(i+1,1)-dens(i,1));beta(1,i)=(-drodh/dens(i,2));
    V(1,i)=Ve*exp(0.01*dens(i,2)*Sw*0.47/(2*beta(1,i)*sin(lambda_entrada)*ro*4*R*R*R*pi/3));
     %V(1,i)=Ve*exp(dens(1,2)*exp(-beta(i)*h(1,i))*Sw*0.47/(2*beta(i)*sin(lambda_entrada)*ro*4*R*R*R*pi/3));
end
h=h(1,1:H/100);h=-h+H;V=fliplr(V);%h=fliplr(h);
plot(h,V),grid on;
Vh=[h;V];clear i;
save Vh
%plot(1:1:T/incr_t,atmosfera_isa)
run EDO_Contorno