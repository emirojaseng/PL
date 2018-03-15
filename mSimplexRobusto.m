function [x0, z0, ban, iter] = mSimplexRobusto(A, b, c)
%purpose: Versión del Simplex más Robusto
% minimizar c^T x
% sujeto a Ax <= b , x >= 0 , b en R^m
%
% In : A ... mxn matrix
% b ... column vector with as many rows as A
% c ... column vector with as many entries as one row of A
%
% Out: xo ... SFB óptima del problema
% zo ... valor óptimo del problema
% ban ... indica casos:
% -1 ... si el conjunto factible es vacio
% 0 ... si se encontro una solución óptima
% 1 ... si la función objectivo no es acotada.
% iter ... es el numero de iteraciones (cambios de variables basicas)
% que hizo el método

if(all(b>=0)) 
    %caso trivial, pues 0 está en el Conjunto Factible
    [x0, z0, ban, iter]= mSimplexFaseII (A,b,c);
    
else
    %caso complicado, hacemos un problema auxiliar
    iter=0;
    sz= size(A);
    m=sz(1);
    n=sz(2);
    c=[c;zeros(m,1)];
    %Checa que no haya renglon A(i,:)>0 para cada b(i)<0
    %falta
    %
    
    %paso 1 variables de holgura
    A=[A, eye(m),-ones(m,1)];
    z=[zeros(n+m,1);1];
    N=[1:n,n+m+1];
    B=n+1:n+m;
    x=zeros(n+m+1,1);
    x(B)=b; %solución basica, no factible
    
    %paso 2 x0 entra, x(n+s) sale
    e=n+m+1; %entra x0
    s=B(find(b==min(b),1)); %sale x(n+s)
    bs=min(b);
    
    %Paso 3: Actualizando a la primera SBF
    x(n+m+1)=-bs; %x0
    x(B)=b-bs;    %xs=0 y el resto de variables se actualizan
    B(B==s)=e;    %se actualizan los indices
    N(N==e)=s;
    
    %Paso 4 Fase II del problema auxiliar
    [x, z0, ban, iter, B, N]=FaseII (A,b,z, x, B, N,n,m, iter);
    
    
    %Checamos los casos
    if(z0==0)
    %x es una SBF para el problema original
            if(~isempty(find(B==n+m+1, 1)))
                %Necesitamos cambiar x0 a las variables no básicas
                %cambiamos variable
            end
    %Buscamos la solución del problema original
    
        N=setdiff(N,[n+m+1]); %Quitamos x0 de los índices
        [x, z0, ban, iter, B, N]=FaseII (A,b,c, x, B, N,n,m, iter);

            %Quitamos las variables de holgura
                    i=1;
                    Baux=B;
                    while i<=length(Baux)
                        if (Baux(i)>n)
                            Baux(i)=[];
                        else
                            i=i+1;
                        end
                    end
                    Baux=sort(Baux);
                    x0=zeros(n,1);
                    x0(Baux)=x(Baux);

    else
    %Caso donde x0>0 y por lo tanto Cf es vacio
            ban=-1;
            x0=NaN;
            z0=NaN;
    end
end
end







function [x, z0, ban, iterAct, B, N]= FaseII (A,b,c, x, B, N,n,m, iter)
%purpose: Fase II empieza con una SBF y encuentra la SBF óptima junto con 
% N y B
% minimizar c^T x
% sujeto a Ax= b , x >= 0 , b en R^m
%
% In : A ... mx(n+m) matrix
% b ... column vector with as many rows as A
% c ... column vector with as many entries as one row of A
% x ... una SBF
% B ... índices de variables básicas
% N ... índices de variables no básicas
% n ... variables en problema original
% m ... restriciones en problema original
%
% Out: xo ... SFB óptima del problema
% zo ... valor óptimo del problema
% ban ... indica casos:
% 0 ... si se encontro una solución óptima
% 1 ... si la función objectivo no es acotada.
% iter ... es el numer´o de iteraciones (cambios de variables basicas)
% que hizo el método

ban=0;
ABInv=inv(A(:,B));
%h conserva el significado que en las notas AB^-1*b
h=ABInv*b;
%r es el vector de Ahorros relativos
r=c(B)'*ABInv*A(:,N)-c(N)';

%Si el vector de ahorros relativos es menor o igual a 0 ya no se entra
%en el ciclo. 
    while (~isempty(find(r>0, 1)))
        
        e=N(find(r==max(r),1));
        He= ABInv*A(:,e);
       %Se comprueba que existe al menos una entrada en He mayor a 0
        if(~isempty(find(He>0,1)))
            aux=h./He;
            aux(aux<=0)=inf; 
            s=B(find(aux==min(aux),1));
            %Acutalizamos los nuevos índices que pertenecen (entrada y
            %salida) de las variables básicas y no básicas
            B(B==s)=e;
            N(N==e)=s;
            ABInv=actInverseR1(A,B,b,x(B),ABInv,A(:,e)-A(:,s),A(:,B(B==e)));
            x(N)=0;
            r=c(B)'*ABInv*A(:,N)-c(N)';
            h=ABInv*b;
            x(B)=h;
            
        else
            ban=1;
            break;
        end
        iter=iter+1;
    end
    
    if(ban==0)
        z0=c(B)'*x(B);
    else
        z0=-inf;
    end
    iterAct=iter;
end

%Función que actualiza la inversa usando Sherman-Morrison
%Reinvierte en caso de ser necesario
function [InvAct]=actInverseR1(A,B,b,h,AInv,u,v)
    tol=0.00001;
    if(norm(b-A(:,B)*h,inf)<tol)
        aux= (AInv*u*v'*AInv)/(1+v'*AInv*u);
        InvAct=AInv-aux;
    else
        InvAct=inv(A(:,B));
    end
end

