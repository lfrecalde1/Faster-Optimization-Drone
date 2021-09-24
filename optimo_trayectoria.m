%% Programa de control predictivo usando casadi como un solver

%% Eliminar variables del sistema
clear all;
close all;
clc;

%% Generacion de los tiempos del sistema
ts = 0.1; 
tfinal = 50;
t = [0:ts:tfinal];

%% Definicion del horizonte de prediccion
N = 50; 

%% Dimension;
a=0.1;

%% Restricciones acciones de control
bounded =[1; -1; pi; -pi];

%% Matrices del sistema
Q=zeros(3,3);
R =zeros(2,2);
Q(1,1) = 5;
Q(2,2) = 5;
Q(3,3) = 0;
R(1,1) = 0.1;
R(2,2) = 0.5;

%% Generacion del soolver por optimizacion
[f, solver, args] = optimization_trayectory(bounded, N, a, ts, Q, R);

%% Definicion de los estados iniciales del sistema
x = -1;
y = -1;
th(1,1) = (180*pi/180);

%% Definicion de la condicion inicial del sistema
hx(1,1)=x+a*cos(th(1));
hy(1,1)=y+a*sin(th(1));

%% Definicion del vector de control inicial
v = zeros(N,2);  % two control inputs 
H0 = repmat([hx(1,1);hy(1,1);th(1,1)],1,N+1)';

%% SENALES DESEADAS DEL SISTEMA
hxd=1*cos(0.5*t);
hxdp=-(0.5)*sin(0.5*t);
hxdpp=-(0.5)*(0.5)*cos(0.5*t);

hyd=1*sin(0.5*t);
hydp=(0.5)*cos(0.5*t);
hydpp=-(0.5)*(0.5)*sin(0.5*t);


hthd=(atan2(hydp,hxdp));

%% GENERACION DE LAS ACCIONES DE CONTROL DE REFERENCIA
vRef = sqrt(hxdp.^2+hydp.^2);
wRef = (hxdp.*hydpp-hydp.*hxdpp)./(hxdp.^2+hydp.^2);

for k=1:length(t)-N

    %% Generacion del vector de estados deseados
    hd=[hxd(1,k);hyd(1,k);hthd(1,k)];
    
    %% deficnion del vector de estados del sistema
    h=[hx(1,k);hy(1,k);th(1,k)];
    
    %% Generacion del; vector de error del sistema
    he(:,k)=hd-h;
    
    args.p(1:3) = h; % initial condition of the robot posture
    for i = 1:N %new - set the reference to track
        %args.p(3*i+1:3*i+3)=[hxd(k+i);hyd(k+i);hthd(k+i)];
        args.p(5*i-1:5*i+1)=[hxd(k+i);hyd(k+i);hthd(k+i)];
        args.p(5*i+2:5*i+3)=[vRef(k+i);wRef(k+i)];
    end 
    args.x0 = [reshape(H0',3*(N+1),1);reshape(v',2*N,1)]; % initial value of the optimization variables
    tic;
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    toc
    sample(k)=toc;
    opti = reshape(full(sol.x(3*(N+1)+1:end))',2,N)';
    H0 = reshape(full(sol.x(1:3*(N+1)))',3,N+1)';
    qpref=[opti(1,1);opti(1,2)];
    %% Resultado del sistema de optimizacion
    u(k) =qpref(1);
    w(k) =qpref(2);
    
    %% Simulacion del sistema
    h=h+system(h,qpref,f,ts);
    
    %% Actualizacion de los estados del sistema
    hx(1,k+1) =h(1);
    hy(1,k+1) =h(2);
    th(1,k+1) =h(3);
    
    %% Actualizacion de los resultados del optimizador para tener una soluciona aproximada a la optima
    v = [opti(2:end,:);opti(end,:)];
    H0 = [H0(2:end,:);H0(end,:)];
end


figure(1)
plot(hx,hy,'b-')
grid on
hold on
plot(hxd,hyd,'g-')

figure(2)
subplot(2,1,1)
plot(t(1:length(u)),u,'b-')
grid on;
hold on;
subplot(2,1,2)
plot(t(1:length(u)),w,'r-')
grid on 
hold on

figure(3)
plot(t(1:length(u)),sample,'b-')
grid on;
hold on;

figure(4)
plot(t(1:length(u)),he(1,:),'b-')
grid on;
hold on;
plot(t(1:length(u)),he(2,:),'r-')