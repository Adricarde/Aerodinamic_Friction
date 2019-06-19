function [qt,h]=calor(h,Vh,dens,M,s)
c=1;
for i=1:size(Vh,2)
    if Vh(1,i)<=h&&c==1;
        V=(Vh(2,i+1)-Vh(2,i))/(Vh(1,i+1)-Vh(1,i))*h+Vh(2,i);
        t=i;c=0;
    end
end
for i=1:M
    if (i-1)*2*pi/M<pi/2||(i-1)*2*pi/M>3*pi/2;
        Re=(i-0.9)*0.1/(2*pi*M)*V*dens(t,2)/(1.9*10^-5);
        Cf=0.664/(Re^0.5);
        qt(i)=0.25*V*V*V*Cf;
    else
        qt(i)=0;
    end
end
h=h-V*s;
end

