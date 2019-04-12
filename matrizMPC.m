function [Psi,Gamma,Tetha,Qex,Rex] =  matrizMPC(A,B,Q,R,S,Hp,Hu)
%% Matriz Psi
Psi=[];
CA=A;
for n= 1:Hp
        Psi=vertcat(Psi,CA);
        CA=A*CA;
end

%% Matriz GAMMA
Ae=[];
Gamma=[];
for i=0:1:Hp-1
    Aex=(A^i)*B;
    Ae=horzcat(Ae,Aex);
    Ai=sum(Ae,2);
    Gamma=vertcat(Gamma,Ai);
end

%% Matriz TETHA
m=size(B,2);
n=length(A);%Dimensi?n de los estados
Tetha=zeros(n*Hp,m*Hu);

for i=1:Hu
    nc=n*(i-1);
    Tetha(nc+1:end,m*(i-1)+1:i*m)=Gamma(1:end-nc,1:m);
end

%% Matriz Q expandida
O=size(S,1);
Aq=eye(Hp);
Qex=kron(Aq,Q);
dQ=size(Qex,1);
Sa=kron(Aq,S);
for i=(dQ-O)+1:size(Qex,1)
    for j=(dQ-O)+1:size(Qex,1)
        Qex(i,j)=Sa(i,j);
    end
end

%% Matriz R expandida
m=size(B,2);
Ar=eye(m*Hu);
Rex=kron(Ar,R);

end

