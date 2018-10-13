clear; clc;
close all;

%% Define symbolic variables
% Parameters
syms K_hsf b_hsf G Jm_R rw Jw m c1 c2 c3 fRLz0 fRRz0

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

% Tire forces 5Burckhardt model)
muRLx = c1*(1-exp(-c2*sRLx)) - c3*sRLx;
muRRx = c1*(1-exp(-c2*sRRx)) - c3*sRRx;

fRLx = muRLx * fRLz;
fRRx = muRRx * fRRz;

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





