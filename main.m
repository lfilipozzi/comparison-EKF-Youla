clear; clc;
close all;

%% Load vehicle parameters
run parameters.m

%% Define initial conditions
U_0 = 100/3.6;
disp(['Initial velocity: ',num2str(3.6*U_0),' km/h'])
chassis_state_init   = [U_0 0 0 U_0/rw U_0/rw U_0/rw U_0/rw 0 0 0];
% rear_axle_state_init = [G*U_0*rw U_0*rw U_0*rw 0 0 0];
rear_axle_state_init = [G*U_0/rw 0];

%% Run simulation
sim('model')