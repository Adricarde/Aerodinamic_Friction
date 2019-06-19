%Calentamiento de una esfera
clear all
close all
load Vh
load Calor
tic
qt=-qt./kterm;
p=a/(N);q=b/(M);s=incr_t;
U=zeros(N,M,L);U(:,:,:)=Tamb;
k=1;
[F,U]=discresp(U,N,M,k,c,p,q,a,qt(:,k));
k=2;  %Iniciamos con Euler
for i=2:N-1;
    for j=2:M-1;
        U(i,j,k)=U(i,j,k-1)+s*F(i-1,j-1);
    end
end

kn=0;
for k=3:5;
    
    for i=N-1:-1:2;
        for j=2:M-1;
            %Adam Bashford 2 para iniciar hasta 5
            [Fn,U]=discresp(U,N,M,k-1,c,p,q,a,qt(:,k-1));
            [Fn0,U]=discresp(U,N,M,k-2,c,p,q,a,qt(:,k-2));
            U(i,j,k)=U(i,j,k-1)+s*0.5*(3*Fn(i-1,j-1)-Fn0(i-1,j-1));
            [Fn1,U]=discresp(U,N,M,k,c,p,q,a,qt(:,k));
            U(i,j,k)=U(i,j,k-1)+s*0.5*(Fn(i-1,j-1)+Fn1(i-1,j-1));
        end
    end
end
for k=5:L;
    
    for i=N-1:-1:2;
        for j=2:M-1;
            %Predictor Adam Bashford
            [Fn3,U]=discresp(U,N,M,k-1,c,p,q,a,qt(:,k-1));
            [Fn2,U]=discresp(U,N,M,k-2,c,p,q,a,qt(:,k-2));
            [Fn1,U]=discresp(U,N,M,k-3,c,p,q,a,qt(:,k-3));
            [Fn0,U]=discresp(U,N,M,k-4,c,p,q,a,qt(:,k-4));
            U(i,j,k)=U(i,j,k-1)+s*(-9*Fn0(i-1,j-1)+37*Fn1(i-1,j-1)-59*Fn2(i-1,j-1)+55*Fn3(i-1,j-1))/24;
            %Corrector Adam Moulton
            [Fn4,U]=discresp(U,N,M,k,c,p,q,a,qt(:,k));
            U(i,j,k)=U(i,j,k-1)+s*(Fn1(i-1,j-1)-5*Fn2(i-1,j-1)+19*Fn3(i-1,j-1)+9*Fn4(i-1,j-1))/24;
        end
    end
    kn=kn+1;
    if kn==10
        disp('Progreso (%)')
        disp(100*k/L)
        kn=0;
    end
end
toc
%Ploteo en polares
figure('units','normalized','outerposition',[0 0 1 1])
rep=1;%v=300:20:400;v=[v 400:100:max(max(max(U)))];
scrsz=get(0,'Screensize');
%figure('Position',[1 scrsz(4)*0.5 scrsz(3) scrsz(4)*0.5])
%figure('Position',[1 1 scrsz(3)*0.4 scrsz(4)*0.4])
while rep==1;
    %figure(1)
    subplot(1,2,2), plot(t,Uh,'-g',t,atmosfera_isa),grid on,xlabel('Tiempo(s)'),ylabel('Altura(m)'),title('Altura en función del tiempo'),legend('Recorrido por hacer','Tropopausa','Estratosfera','Estratosfera','Estratopausa','Mesosfera','Mesosfera','Mesopausa');
    
    for k=1:L;
        for i=1:N;
            for j=1:M;
                x(i,j)=i*p*cos((j)*q);
                y(i,j)=i*p*sin((j)*q);
                T(i,j)=U(i,j,k);
            end
        end
        for i=1:N;
            for j=M:2*M;
                x(i,j)=i*p*cos((j)*q);
                y(i,j)=i*p*sin((j)*q);
            end
        end
        T=[T, fliplr(T)];
        
        %figure(1)
        subplot(1,2,1), surf(x,y,T),axis([-a a -a a 300 max(max(max(U)))]),view(0,90);colormap hot
        subplot(1,2,2),
        hold on
        plot(t(k),Uh(k),':ok');
        hold off
        pause(0.01);
        
%         figure(2)
%         h=pcolor(x,y,T);
%         set(h,'EdgeColor','none');colormap hsv;colorbar
%         hold on
%         contour(x,y,T,v,'b','ShowText','on')
%         hold off
        
        clear T
    end
    rep=input('Escriba 1 si quiere repetir el plot, si no, escriba 0 \n');
end

clc
