function [Fv,qt]=fvel(U,i,Ve,H)
load Vh
dens=flipud(dens);c=1;
%R=0.1;Sw=2*pi*R*R;%Metros, Metros/Segundo, Metros, Metros*Metros, 
%ro=2700;lambda_entrada=-pi/2;%Kilos/metro,radianes
%mu=1.9*10^-5;%viscosidad dinámica
for j=1:size(Vh,2);
    if Vh(1,j)<=U(1,i)&&c==1;
        c=0;t=j;
    end
end
dens=dens(1:H/100+1,1:2);
drodh=(dens(t+1,2)-dens(t,2))/(dens(t+1,1)-dens(t,1));beta=-drodh/dens(t,2);
Fv=-Ve*exp(dens(t,2)*Sw*0.47/(2*beta*sin(lambda_entrada)*ro*4*R*R*R*pi/3));
for j=1:M
    if (j-1)*2*pi/(M-1)<pi/2
        
        Re=(j-0.8)*R/(2*pi*M)*Fv*dens(t,2)/(1.9*10^-5);
        Cf=0.664/(abs(Re)^0.5);
        qt(j)=-0.25*Fv*Fv*Fv*Cf*dens(t,2);
    else
        qt(j)=500.5;
    end
end
end
