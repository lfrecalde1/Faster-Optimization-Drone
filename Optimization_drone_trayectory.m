%% Programa de control predictivo para un drone basado en optimizacion usando Casadi

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

%% Definicion de las constantes del sistema
a = 0.0;
b = 0.0;
L =[a;b];

%% Definicion de los limites de las acciondes de control
bounded = [1.2; -1.2; 1.2; -1.2; 1.2; -1.2; 1.5; -1.5];

%% Seccion para cargar los parametros dinamicos del sistema
load('DINAMICA_DRONE.mat')
chi_real=x'*ones(1,length(t));

%% Deficion de la matriz de la matriz de control
Q = 0*eye(4);
Q(1,1) = 1;
Q(2,2) = 1;
Q(3,3) = 1;
Q(4,4) = 0.1;

%% Definicion de la matriz de las acciones de control
R = 0.2*eye(4);


%% Definicion de los estados iniciales del sistema
x = 0;
y = 0;
z = 0;
th(1,1) = (180*pi/180);

%% Defincion de la condicion inial del sistema
hx(1,1) = x+a*cos(th(1))-b*sin(th(1));
hy(1,1) = y+a*sin(th(1))+b*cos(th(1));
hz(1,1) = z;

%% Definicion del vectro de control inicial del sistema
v = zeros(N,4);
H0 = repmat([hx(1,1);hy(1,1);hz(1,1);th(1,1)],1,N+1)';

%% DEFICNION DE LAS VELOCIDADES DE CONTROL REALES 
ul(1)=0;
um(1)=0;
un(1)=0;
w(1)=0;

%% Senales deseadas del sistema

%% Senals deseadas en x
hxd=0.25*t+2;
hxdp=0.25*ones(1,length(t));
hxdpp=0*ones(1,length(t));

hyd=2*sin(t/8)+0.05*t-4;
hydp=(1/8)*2*cos(t/8)+0.05;
hydpp=-(1/8)*(1/8)*2*sin(t/8);

hzd=10+1.5*sin(t/10);
hzdp=1.5*(1/10)*cos(t/10);
hzdpp=-1.5*(1/10)*1/10*sin(t/10);

hthd=(atan2(hydp,hxdp));
hthdp=diff([0 hthd])/ts;

%% Seccion para cargar los parametros dinamicos del sistema
load('DINAMICA_DRONE.mat')
chi_real=x'*ones(1,length(t));

%% Obstacle del sistema
ubicacion= 250
obs = [hxd(ubicacion);hyd(ubicacion);hzd(ubicacion)];
 
%% Definicion del optimizador
[f, solver, args] = mpc_drone(bounded, N, L, ts, Q, R, obs);

% Simulacion del sistema
for k=1:length(t)-N

    %% Generacion del vector de estados deseados
    hd=[hxd(1,k);hyd(1,k);hzd(1,k);hthd(1,k)];
    
    %% deficnion del vector de estados del sistema
    h=[hx(1,k);hy(1,k);hz(1,k);th(1,k)];  
    
    %% Generacion del; vector de error del sistema
    he(:,k)=hd-h;
    
    %% Generacion del vector de velocidades reales del sistema
    v_real=[ul(k);um(k);un(k);w(k)];
    
    args.p(1:4) = h; % Generacion del estado del sistema
    
    for i = 1:N % z
        args.p(4*i+1:4*i+4)=[hxd(k+i);hyd(k+i);hzd(k+i);hthd(k+i)];
%         args.p(4*i+5:4*i+7)=obs;
    end 
    
    args.x0 = [reshape(H0',4*(N+1),1);reshape(v',4*N,1)]; % initial value of the optimization variables
    tic;
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    toc
    sample(k)=toc;
    opti = reshape(full(sol.x(4*(N+1)+1:end))',4,N)';
    H0 = reshape(full(sol.x(1:4*(N+1)))',4,N+1)';
    qpref=[opti(1,1);opti(1,2);opti(1,3);opti(1,4)];
    
    %% Resultado del sistema de optimizacion
    ul_c(k) = qpref(1);
    um_c(k) = qpref(2);
    un_c(k) = qpref(3);
    w_c(k) = qpref(4);
    
    %% Dinamica del sistema 
    DINAMICA_DRONE = Dinamica(qpref,v_real,chi_real(:,k),ts);
    ul(k+1)=DINAMICA_DRONE(1);
    um(k+1)=DINAMICA_DRONE(2);
    un(k+1)=DINAMICA_DRONE(3);
    w(k+1)= DINAMICA_DRONE(4);
    
    %% Simulacion del sistema
    h=h+system(h,[ul(k+1);um(k+1);un(k+1);w(k+1)],f,ts);
    
    %% Actualizacion de los estados del sistema
    hx(1,k+1) = h(1);
    hy(1,k+1) = h(2);
    hz(1,k+1) = h(3);
    th(1,k+1) = h(4);
    
    %% Actualizacion de los resultados del optimizador para tener una soluciona aproximada a la optima
    v = [opti(2:end,:);opti(end,:)];
    H0 = [H0(2:end,:);H0(end,:)];
end
close all; paso=1; 
%a) Parámetros del cuadro de animación
figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 8 3]);
h = light;
h.Color=[0.65,0.65,0.65];
h.Style = 'infinite';
%b) Dimenciones del Robot
   Drone_Parameters(0.02);
%c) Dibujo del Robot    
    G2=Drone_Plot_3D(hx(1),hy(1),hz(1),0,0,th(1));hold on

    plot3(hx(1),hy(1),hz(1),'--','Color',[56,171,217]/255,'linewidth',1.5);hold on,grid on   
    plot3(hxd(1),hyd(1),hzd(1),'Color',[32,185,29]/255,'linewidth',1.5);
    plot3(hxd(ubicacion),hyd(ubicacion),hzd(ubicacion),'*r','linewidth',1.5);
view(20,15);

for k = 1:10:length(t)-N
    drawnow
    delete(G2);
   
    G2=Drone_Plot_3D(hx(k),hy(k),hz(k),0,0,th(k));hold on
    
    plot3(hxd(1:k),hyd(1:k),hzd(1:k),'Color',[32,185,29]/255,'linewidth',1.5);
    plot3(hx(1:k),hy(1:k),hz(1:k),'--','Color',[56,171,217]/255,'linewidth',1.5);
    
    legend({'$\mathbf{h}$','$\mathbf{h}_{des}$'},'Interpreter','latex','FontSize',11,'Location','northwest','Orientation','horizontal');
    legend('boxoff')
    title('$\textrm{Movement Executed by the Aerial Robot}$','Interpreter','latex','FontSize',11);
    xlabel('$\textrm{X}[m]$','Interpreter','latex','FontSize',9); ylabel('$\textrm{Y}[m]$','Interpreter','latex','FontSize',9);zlabel('$\textrm{Z}[m]$','Interpreter','latex','FontSize',9);
    
end

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);

plot(t(1:length(he)),he(1,:),'Color',[226,76,44]/255,'linewidth',1); hold on;
plot(t(1:length(he)),he(2,:),'Color',[46,188,89]/255,'linewidth',1); hold on;
plot(t(1:length(he)),he(3,:),'Color',[26,115,160]/255,'linewidth',1);hold on;
plot(t(1:length(he)),he(4,:),'Color',[83,57,217]/255,'linewidth',1);hold on;
grid on;
legend({'$\tilde{h_{x}}$','$\tilde{h_{y}}$','$\tilde{h_{z}}$','$\tilde{h_{\psi}}$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Evolution of Control Errors}$','Interpreter','latex','FontSize',9);
ylabel('$[m]$','Interpreter','latex','FontSize',9);


figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1:length(ul_c)),ul_c,'Color',[226,76,44]/255,'linewidth',1); hold on
plot(t(1:length(ul_c)),um_c,'Color',[46,188,89]/255,'linewidth',1); hold on
plot(t(1:length(ul_c)),un_c,'Color',[26,115,160]/255,'linewidth',1); hold on
plot(t(1:length(ul_c)),w_c,'Color',[83,57,217]/255,'linewidth',1); hold on
plot(t(1:length(ul)),ul,'--','Color',[226,76,44]/255,'linewidth',1); hold on
plot(t(1:length(ul)),um,'--','Color',[46,188,89]/255,'linewidth',1); hold on
plot(t(1:length(ul)),un,'--','Color',[26,115,160]/255,'linewidth',1); hold on
plot(t(1:length(ul)),w,'--','Color',[83,57,217]/255,'linewidth',1); hold on
grid on;
legend({'$\mu_{lc}$','$\mu_{mc}$','$\mu_{nc}$','$\omega_{c}$','$\mu_{l}$','$\mu_{m}$','$\mu_{n}$','$\omega$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Control Values}$','Interpreter','latex','FontSize',9);
ylabel('$[rad/s]$','Interpreter','latex','FontSize',9);
xlabel('$\textrm{Time}[s]$','Interpreter','latex','FontSize',9);

figure
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [4 2]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 10 4]);
plot(t(1:length(sample)),sample,'Color',[46,188,89]/255,'linewidth',1); hold on
grid on;
legend({'$t_{sample}$'},'Interpreter','latex','FontSize',11,'Orientation','horizontal');
legend('boxoff')
title('$\textrm{Sample Time}$','Interpreter','latex','FontSize',9);
ylabel('$[s]$','Interpreter','latex','FontSize',9);
xlabel('$\textrm{Time}[s]$','Interpreter','latex','FontSize',9);