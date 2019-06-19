
function [FU,U]=discresp(U,N,M,k,c,p,q,a,qt)
%p es el incremento de r y q es el incremento de fi
%U(N,fi_ini:fi_fin,k)=Tamb+frozamiento; %Zona afectada por la fricción del aire
r1=c*c/(p*p); %simplifiación del coeficiente de calor para cada material con el incremento de r
for i=N:-1:1;
    for j=1:M;
      if i==N
         U(i,j,k)=(-qt(j)*2*p+4*U(i-1,j,k)-U(i-2,j,k))/3;      
      elseif i==1;
         U(i,j,k)=(4*U(i+1,j,k)-U(i+2,j,k))/3;
     end
      if j==1 && i<N;       
       U(i,1,k)=(4*U(i,2,k)-U(i,3,k))/3;
       U(i,M,k)=(4*U(i,M-1,k)-U(i,M-2,k))/3;
      end
    end
end
            

for i=N-1:-1:2;
    for j=2:M-1;
        FU(i,j)=(r1/(i*i))*(i*(U(i+1,j,k)-U(i-1,j,k))+i*i*(U(i+1,j,k)+U(i-1,j,k)-2*U(i,j,k))+((U(i,j+1,k)+U(i,j-1,k)-2*U(i,j,k))/(q*q)));
    end
end
FU=FU(2:N-1,2:M-1);
end
%La ecuación del calor en esféricas la he desarrolado un poco para que
%quedara más sencilla, pero más larga.

