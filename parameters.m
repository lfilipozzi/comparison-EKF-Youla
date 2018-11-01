%% General parameters
g = 9.81;

%% Chassis parameters
m  = 1370+2*80; % Vehicle mass (kg)
Jy = 2315.3;    % Yaw moment of inertia (kg/m^2)
l  = 2.78;      % Wheelbase (m)
a  = 1.1;       % Distance from the front axle to the CG (m)
b  = l - a;     % Distance from the rear axle to the CG (m)
T  = 1.55;      % Track width (m)
h  = 0.520;     % Height of the CG (m)

% Compute static normal forces and weight transfer due to acceleration
fFLz0 = m*g*b/(2*l); % Static normal force on the FR wheel (N)
fFRz0 = m*g*b/(2*l); % Static normal force on the FL wheel (N)
fRLz0 = m*g*a/(2*l); % Static normal force on the RL wheel (N)
fRRz0 = m*g*a/(2*l); % Static normal force on the RR wheel (N)
DFx   = m*h/(2*l);   % Level arm due to longitudinal acceleration (m)
DFFy  = m*h*b/(T*l); % Level arm due to lateral acceleration front axle (m)
DFRy  = m*h*a/(T*l); % Level arm due to lateral acceleration rear axle (m)

%% Wheels parameters
Jw = 0.9;       % Wheel inertia (kg.m^2)
rw = 0.325;     % Wheel radius (m)

%% Rear axle parameters
% Output gearing (Ford parameters)
N1     = 26;
N_2e   = 83;
N_2m   = 101;%93; % Modified to obtain g_gearbox = 10
N3     = 59; 
N4     = 23;
N5     = 59;
rho_PS = 0.3953;
Gfd_ps = N5/N4;     % Final drive ratio
Gm_ps  = N_2m/N1;   % Power split motor gear ratio
G1_PS  = 2/Gfd_ps; 
G2_PS  = Gm_ps;     % Power split motor gear ratio

% Calculate equivalent driveshaft inertias
J_df  = 0.041937;   % Rear differential housing inertia
J_gen = 0.0106;     % PS generator inertia
J_eng = 0.18;       % PS engine inertia
J_mf  = 0.0351;     % PS motor inertia
J_eng_eq =  J_eng + J_gen*((1+rho_PS)/rho_PS)^2;

% Lump inertias of PS motor and front diff into layshaft inertia
m_lay = 5;          % Equivalent gear masses (GUESS)
r_lay = 100/1000;   % Equivalent gear radius (GUESS)
J_lay = 0.5*m_lay*r_lay^2; % Inertia of PS layshaft (GUESS)
J_ps = J_lay + J_df*Gfd_ps^2 + ...    % Equivalent inertia of PS layshaft
    J_mf*G2_PS^2 + J_eng_eq*(N_2e/N1)^2; 

J_ps = J_ps/10;

G = Gm_ps * Gfd_ps; % Global ratio
Jm_R = J_mf;

Gring = Gm_ps;

%% Lumped halfshaft
% Rear axle
K_hsf_RL = 80*180/pi;   % RL halfshaft compliance (Nm/rad)
b_hsf_RL = 0.8*180/pi;  % RL halfshaft damping (Nms/rad)
K_hsf_RR = K_hsf_RL;    % RR halfshaft compliance (Nm/rad)
b_hsf_RR = b_hsf_RL;    % RR halfshaft damping (Nms/rad)

% Combined stiffness and damping for RL and RR in parallel
K_hsf = 1/(1/K_hsf_RL + 1/K_hsf_RR);
b_hsf = 1/(1/b_hsf_RL + 1/b_hsf_RR);

J_hsf = 0.009;          % Halfshaft inertia

% Rear axle motor
K_sf_R = 200*180/pi;   % Motor axle halfshaft compliance (Nm/rad)
b_sf_R = 0.8*180/pi;   % Motor axle halfshaft damping  (Nms/rad)

%% Tire parameters
% Dugoff tire model
Cx = 866.05/0.02;           % Tire longitudinal stiffness
Cy = 248.08/(0.5*pi/180);   % Tire lateral stiffness
epsDugoff = 0.01;

% Burckhardt tire model
c1 = 1.2801;
c2 = 23.99;
c3 = 0.25;

% Pacejka's magic formula
B = 7.1;
C = 1.3;
D = 1;

% % Coef from Mathworks:
% % https://www.mathworks.com/help/physmod/sdl/ref/tireroadinteractionmagicformula.html
% % Dry tarmac	
% B = 10;	
% C = 1.9;
% D = 1;
% E = 0.97;
% % Wet tarmac
% B = 12;
% C = 2.3;
% D = 0.82;
% E = 1;
% % Snow	
% B = 5;
% C = 2;
% D = 0.3;
% E = 1;
% % Ice
% B = 4;
% C = 2;
% D = 0.1;
% E = 1;

%% Road parameters
alpha = 0;  % Road inclination (%)
alpha = atan(alpha/100);    % Convert into angle

