clear; clc;
close all;

%% Define symbolic variables
% Parameters
% Do not use <variable>_est here only <variable>
syms K_hsf b_hsf G Jm_R rw Jw m c1 c2 c3 fRLz0 fRRz0 ...
    muRL0 muRR0 CRx CRy epsDugoff

% States
syms theta_hsf wm_R ww_RL ww_RR U

% Inputs
syms tau_m_R

%% Equation of motion
% Tire slip
sRLx = (rw*ww_RL - U) / U;
sRRx = (rw*ww_RR - U) / U;

% % Normal forces (TODO: solve algebraic loop)
fRLz = fRLz0;% + DFx*ax;
fRRz = fRRz0;% + DFx*ax;

% % Tire forces (Burckhardt model)
% muRLx = c1*(1-exp(-c2*sRLx)) - c3*sRLx;
% muRRx = c1*(1-exp(-c2*sRRx)) - c3*sRRx;
% fRLx = muRLx * fRLz;
% fRRx = muRRx * fRRz;


% Tire forces (Dugoff)
alphaRL = 0;
alphaRR = 0;
% [fRLx,fRLy] = Dugoff(CRx,CRy,fRLz,sRLx,alphaRL,U,muRL0,epsDugoff);
% [fRRx,fRRy] = Dugoff(CRx,CRy,fRRz,sRRx,alphaRR,U,muRR0,epsDugoff);
% Dugoff function use if condition (ccanot be used with symbolic variable)

% % Consider design in the linear region (lambda > 1 in Dugoff.m file)
% disp('Design in the linear region (uncomment lines to obtain matrix in non-linear region)')
% fRLx = CRx * sRLx / abs(1-abs(sRLx));
% fRRx = CRx * sRRx / abs(1-abs(sRRx));

% % Consider design in the non-linear region (lambda <= 1 in Dugoff.m file)
% disp('Design in the non-linear region for RL (uncomment lines to obtain matrix in linear region)')
% % RL
% fRLxd = CRx * sRLx / abs(1-abs(sRLx));
% fRLyd = 0;
% lambdaRL = muRL0 * (1 - epsDugoff*U*abs(sRLx)) / ...
%     (2 * sqrt((fRLxd/fRLz)^2 + (fRLyd/fRLz)^2));
% fRLx = fRLxd * 2 * lambdaRL * (1 - lambdaRL/2);
% % RR
% fRRx = CRx * sRRx / abs(1-abs(sRRx));

% % Consider design in the non-linear region (lambda <= 1 in Dugoff.m file)
% disp('Design in the non-linear region for RR (uncomment lines to obtain matrix in linear region)')
% % RL
% fRLx = CRx * sRLx / abs(1-abs(sRLx));
% % RR
% fRRxd = CRx * sRRx / abs(1-abs(sRRx));
% fRRyd = 0;
% lambdaRR = muRR0 * (1 - epsDugoff*U*abs(sRRx)) / ...
%     (2 * sqrt((fRRxd/fRRz)^2 + (fRRyd/fRRz)^2));
% fRRx = fRRxd * 2 * lambdaRR * (1 - lambdaRR/2);

% Consider design in the non-linear region (lambda <= 1 in Dugoff.m file)
disp('Design in the non-linear region (uncomment lines to obtain matrix in linear region)')
% RL
fRLxd = CRx * sRLx / abs(1-abs(sRLx));
fRLyd = 0;
lambdaRL = muRL0 * (1 - epsDugoff*U*abs(sRLx)) / ...
    (2 * sqrt((fRLxd/fRLz)^2 + (fRLyd/fRLz)^2));
fRLx = fRLxd * 2 * lambdaRL * (1 - lambdaRL/2);
% RR
fRRxd = CRx * sRRx / abs(1-abs(sRRx));
fRRyd = 0;
lambdaRR = muRR0 * (1 - epsDugoff*U*abs(sRRx)) / ...
    (2 * sqrt((fRRxd/fRRz)^2 + (fRRyd/fRRz)^2));
fRRx = fRRxd * 2 * lambdaRR * (1 - lambdaRR/2);


% Halshaft torque
tau_hsf = K_hsf * theta_hsf + b_hsf * (2/G*wm_R - ww_RL - ww_RR); 

% Equation of motion
wm_R_dot      = 1/Jm_R * (tau_m_R - 2/G * tau_hsf);  
theta_hsf_dot = 2/G * wm_R - ww_RL - ww_RR;
ww_RL_dot = (tau_hsf - rw * fRLx) /Jw;
ww_RR_dot = (tau_hsf - rw * fRRx) /Jw;
U_dot     = 1/m * (fRLx + fRRx);

% Return states
x = [theta_hsf wm_R ww_RL ww_RR U];

% Return state derivatives
x_dot = [theta_hsf_dot wm_R_dot ww_RL_dot ww_RR_dot U_dot];

% Return inputs
u = tau_m_R;

% Return outputs
y = [U_dot ww_RL ww_RR];

%% Compute Jacobian
A_sym = jacobian(x_dot,x);
B_sym = jacobian(x_dot,u);
C_sym = jacobian(y,x);
D_sym = jacobian(y,u);

A_sym = simplify(A_sym);
B_sym = simplify(B_sym);
C_sym = simplify(C_sym);
D_sym = simplify(D_sym);

%% Evaluate matrices
run parameters.m

U = 50/3.6;
theta_hsf = 1*pi/180;
wm_R = 1*U/rw*G;
ww_RL = 1*U/rw;
ww_RR = ww_RL;

A = eval(A_sym);
B = eval(B_sym);
C = eval(C_sym);
D = eval(D_sym);

rank(obsv(A,C))





