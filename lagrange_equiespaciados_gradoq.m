function [I,xb,fr,Rn,x]=lagrange_equiespaciados_gradoq(N,incr,q,Vh)
c=0;n=N+1;xb(1,1+2/incr)=0;I(1,1+2/incr)=0;x(1,N+2)=0;prod(1,N+2)=1;b(1,N+1)=0;Rn(1,1+2/incr)=0;fr(1,1+2/incr)=0;
ini=1;fin=q+1;
for i=1:N+2;
        x(1,i)=-1+2*(i-1)/(n);
end
for k=-1:incr:1;
c=c+1;

xb(1,c)=k;I(1,c)=0;fr(1,c)=f(xb(1,c));
if xb(1,c)>x(1,fin);
    ini=fin;
    fin=fin+q;
    if fin>(N+2);
        fin=N+2;
    end
end
for j=ini:fin;
    x(1,j)=-1+2*(j-1)/(n);
    b(1,j)=f(x(1,j));
    prod(1,j)=1;
    for i=ini:fin;
        x(1,i)=-1+2*(i-1)/(n);
        if i~=j;
        prod(1,j)=(xb(1,c)-x(1,i))/(x(1,j)-x(1,i))*prod(1,j);
        end
    end
    I(1,c)=I(1,c)+b(1,j)*prod(1,j);
end
Rn(1,c)=abs(fr(1,c)-I(1,c));
end
end