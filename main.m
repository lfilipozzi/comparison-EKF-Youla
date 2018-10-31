clear; clc;
close all;

%% Load vehicle parameters
run parameters.m

Gc = load('Youla_controller.mat','Gc');
Gc = Gc.Gc;

%% Define initial conditions
U_0 = 50/3.6;
disp(['Initial velocity: ',num2str(3.6*U_0),' km/h'])
disp(['Road inclination: ',num2str(100*tan(alpha)),'m / 100m'])
chassis_state_init   = [U_0 0 0 U_0/rw U_0/rw U_0/rw U_0/rw 0 0 0];
rear_axle_state_init = [G*U_0/rw 0];

%% Set sensor parameters
Ts_sensor = 0.001;  % Sampling time of the sensors
sigma_ax  = 0.5;    % Covariance of accelerometers
sigma_wij = 0.1;    % Covariance of wheel speed sensors

%% Set parameters of EKF and UKF
Ts_EKF = 0.001;
Ts_UKF = 0.001;
Qk_EKF = eye(1)/1000000;
Qk_UKF = eye(5)/1000000;
Rk = diag([sigma_ax sigma_wij sigma_wij]);

set_param('model/EKF',...
    'x_init',strcat('[',num2str(0),';',...
        num2str(U_0/rw*G),';',...
        num2str(U_0/rw),';',...
        num2str(U_0/rw),';',...
        num2str(U_0),']'),...
    'Ts'   ,num2str(Ts_EKF),...
    'K_hsf',num2str(K_hsf),...
    'b_hsf',num2str(b_hsf),...
    'G'    ,num2str(G),...
    'Jm_R' ,num2str(Jm_R),...
    'rw'   ,num2str(rw),...
    'Jw'   ,num2str(Jw),...
    'm'    ,num2str(m),...
    'c1'   ,num2str(c1),...
    'c2'   ,num2str(c2),...
    'c3'   ,num2str(c3),...
    'fRLz0',num2str(fRLz0),...
    'fRRz0',num2str(fRRz0))

set_param('model/UKF',...
    'x_init',strcat('[',num2str(0),';',...
        num2str(U_0/rw*G),';',...
        num2str(U_0/rw),';',...
        num2str(U_0/rw),';',...
        num2str(U_0),']'),...
    'Ts'   ,num2str(Ts_UKF),...
    'K_hsf',num2str(K_hsf),...
    'b_hsf',num2str(b_hsf),...
    'G'    ,num2str(G),...
    'Jm_R' ,num2str(Jm_R),...
    'rw'   ,num2str(rw),...
    'Jw'   ,num2str(Jw),...
    'm'    ,num2str(m),...
    'c1'   ,num2str(c1),...
    'c2'   ,num2str(c2),...
    'c3'   ,num2str(c3),...
    'fRLz0',num2str(fRLz0),...
    'fRRz0',num2str(fRRz0))

%% Run simulation
m = 1.2*m;
sim('model')

%% Plot results
figure(1)
hold on
box on
plotErrorPDF(EKF.fRLx,fijx_true.fRLx,1000)
plotErrorPDF(UKF.fRLx,fijx_true.fRLx,1000)
plotErrorPDF(Youla.fRLx,fijx_true.fRLx,1000)
legend('EKF','UKF','Youla')

