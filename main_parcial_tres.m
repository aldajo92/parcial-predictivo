%MPC por Zonas
clc; clear all; close all;

Ts = 0.1;
Tf = 30;

% Systema continuo
A=[0 1;0 0];
B=[0;1];
C=[1 0];
sysC=ss(A,B,C,0);

% Systema discreto
sysD=c2d(sysC,Ts);
Ad=sysD.A;
Bd=sysD.B;       

nes=size(A,1);        % cantidad de estados
men=size(B,2);        % cantidad de entradas
p=size(C,1);          % Salidas controladas

%% Parametros de Control
Hp=20;
Hu=5;
Q=100*eye(p);
R=.0001;
S=1*Q;
X0=[1;-1];             % estado inicial
U=0;                   % u(k-1) ultimo control aplicado
y_Ref = 2;

%% Referencia
Ref = [];
for i = 1:Hp
    Ref = vertcat(Ref,y_Ref);
end

% Matriz Cz (C expandida)
Cz = C_expand(C,Hp);

% Matriz T
IT=eye(p);
T=[];
for i = 1:Hp
    T=vertcat(IT,T);
end

% Matriz Tetha 1
[Psi,Gamma,Tetha,Qex,Rex] =  matrizMPC(Ad,Bd,Q,R,S,Hp,Hu);
Tetha1a=Cz*Tetha;
Tetha1=horzcat(Tetha1a,-T);

% Matriz Rex1
Rex1= [Rex zeros(men*Hu,p);zeros(p,men*Hu) zeros(p)];

% Matriz Hz y Gz
Hz=Tetha1'*Qex*Tetha1 + Rex1;
Gz=2*Tetha1'*Qex;

%% Delta U1

dumin = -1.5;
dumax = 1.5;

dymin=1.8;
dymax=2.3;

% => deltaU
Im = [];
for i = 1:Hu
    Im = vertcat(Im,eye(p));     
end

% => deltaY
Iyd = [];
for i = 1:Hp       
    Iyd = vertcat(Iyd,eye(p));      
end

LiU=[Im*dumin;dymin]; %Limite Inferior para deltaU1
LiS=[Im*dumax;dymax]; %Limite Superior para deltaU1

%% Restricciones para Y

Ymax=[2;5];          %Limite superior para y1, y2
Ymin=[-2;-5];        %Limite inferior para y1, y2

Iy = [];
for i = 1:Hp
    Iy = vertcat(Iy,eye(nes));  
end

FF2=tril(ones(Hp));           % Se encuentra la matriz triangualar inferior
F2=kron(FF2,eye(nes));        % coloca FF2 en identidad tamano n
F3=F2*Tetha;                  % I*Tetha
F31=[F3 zeros(nes*Hp,p)];     % Restriccion para deltaU1

%% Restricciones para U

Umin=-10;
Umax=10;

FF=tril(ones(Hu));            %Se encuentra una matriz triangular inferior
F1=kron(FF,eye(men));         
F11=[F1 zeros(men*Hu,p)];     %Adecuamos la restriccion para deltaU1

%% matriz de limites de U
Im = [];
for i = 1:Hu
    Im = vertcat(Im,eye(men));  
end

%% tiempos de muestreo
time = 0:Ts:Tf;
Samples = size(time, 2);

%% Simulacion
x(:,1)=[1;-1];              %Se define el vector de estados
for j=1:Samples
    E = Ref -(Cz*Psi*x(:,j))-(Cz*Gamma*U);     % Error

    lb = Im*(-Umin+U);  %U inferior
    ub = Im*(Umax-U);   %U superior

    ly= ((Iy+Psi)*(-Ymin+x(:,1))) + (Gamma*U);  %Limite inferior para Y
    uy= ((Iy+Psi)*(Ymax-x(:,1))) + (Gamma*U);   %Limite superior para Y
    
    DeltaUop=quadprog(2*Hz,-Gz*E,[F11;-F11;-F31;F31],[lb;ub;ly;uy],[],[],LiU,LiS);
    
    du(j)=DeltaUop(1);
    
    uv(j)=du(j) + U;   % Accion de control U
    
    x(:,j+1)=Ad*x(:,j)+ Bd*uv(j);
    y(j)=C*x(:,j+1);
    U=uv(j);
end


subplot(2,2,1);
stairs(x(1,:)');grid on;xlim([0 Samples]);title('X_1');

subplot(2,2,2);
stairs(uv);grid on;xlim([0 Samples]);title('Control U');

subplot(2,2,3);
stairs(du);grid on;xlim([0 Samples]);title('Delta U');

subplot(2,2,4);
stairs(y');grid on;xlim([0 Samples]);title('Salida');

