clc
clear all
close all
load Vh
load Calor
N=15;M=15;
p=a/(N);q=b/(M);s=incr_t;

%Error local de truncación
r=0.5;fi=pi/4;ind_p=0;ind_q=0;
for i=0.4:0.05:2
    ind_p=ind_p+1;
    p1=10^(-i);pvect(ind_p)=p1;
    for j=0.4:0.05:2
     ind_q=ind_q+1;
     q1=10^(-j);qvect(ind_q)=q1;
    T(ind_p,ind_q)=c*c*(2*exp(r)*sin(fi)/r+exp(r)*sin(fi)-exp(r)*sin(fi)/(r*r))-(c*c)*((exp(r+p1)*sin(fi)-exp(r-p1)*sin(fi))/r+(exp(r-p1)*sin(fi)+exp(r+p1)*sin(fi)-2*exp(r)*sin(fi))+((exp(r)*sin(fi+q1)+exp(r)*sin(fi-q1)-2*exp(r)*sin(fi))/(q1*q1)) );
    %disp(T)
    T(ind_p,ind_q)=log(T(ind_p,ind_q));
    end
    ind_q=0;
    i/3
end
pvect=fliplr(pvect);qvect=fliplr(qvect);T=T';
%T=log(T);
figure(3)
%plot(pvect,T(:,1))
h=surf(pvect,qvect,T),view(0,90),xlabel('Incremento Radial'),ylabel('Incremento Angular'),title('Error Local (log)');
set(h,'EdgeColor','none');colormap jet;colorbar;
%hold on
%contour(pvect,qvect,T,15,'b','ShowText','on')
%hold off



%Estabilidad
U(N,M,L)=0;
r1=c*c/(p*p); 
A(N-2,M-2)=0;
ind_s=0;
for s=1.5:-0.01:0.01
        U=eye(N);ind_s=ind_s+1;A=0;
for i=2:1:N-1;
    for j=2:M-1;
        FU(i-1,j-1)=(r1/(i*i))*(i*(U(i+1,j,1)-U(i-1,j,1))+i*i*(U(i+1,j,1)+U(i-1,j,1)-2*U(i,j,1))+((U(i,j+1,1)+U(i,j-1,1)-2*U(i,j,1))/(q*q)) );
    end
end
%FU=FU(2:N-1,2:M-1);
A=FU+A;

A(max(N,M)-2,max(N,M)-2)=0;
%A=A(2:N-1,2:N-1);
A=A';
B=eye(size(A,1))+s*A;

svect(ind_s)=log(s);
vectA(ind_s,:)=s*eig(A)';n=1;maxautovalor(ind_s)=max(s*eig(A));
end
%vectstr=num2str(vectA);




figure(1)
%plot(real(vectA),imag(vectA),' ob'),grid on, axis([-1 1 -1 1])%,text(-1, -0.5, vectstr);
plot(svect,maxautovalor),grid on, axis tight;
if max(real(vectA))>0.0001
    disp('Hay autovalores positivos en A')
elseif min(real(vectA))<-1.2
    disp('Hay autovalores erroneos')
else
    disp('No hay autovalores positivos en A')
end
vectB=eig(B)';
vectstr=num2str(vectB);
figure(2)
plot(real(vectB),imag(vectB),' og'),grid on, axis equal,text(0, 0, vectstr);
if max(real(vectB))>1.00001%&&max(vectB)~1;
    disp('Hay autovalores erroneos en B')
else
    disp('No hay autovalores erroneos en B')
end



