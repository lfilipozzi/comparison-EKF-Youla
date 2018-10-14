% Youla Controller Output Observer for 
% Longitudinal Friction estimation


%% Linearizing the Plant 
% Load Parameters
parameters
Jm = Jm_R;
K_hsr = K_hsf;
b_hsr = b_hsf;

% Symbolic 
syms  wm the_hsr wL wR real % States
syms FL FR tau real% Inputs

%% Equations for the Plant
tau_hsr = K_hsr*the_hsr + b_hsr*(2/G * wm - wL - wR);

dwm = (tau-2*tau_hsr/G)/Jm;
dthe_hsr = 2*wm/G - wL - wR;
dwL = (tau_hsr-rw*FL);
dwR = (tau_hsr-rw*FR);

%% States/Inputs/Outputs
dx = [dwm,dthe_hsr,dwL,dwR]';
x = [wm,the_hsr,wL,wR]';
u = [FL,FR]';
y = [wL,wR]';

%% Matrixification
A = double(jacobian(dx,x)); 
B = double(jacobian(dx,u));
C = double(jacobian(y,x));
D = double(jacobian(y,u));

s = tf('s');
G = minreal(C*inv(s*eye(4)-A)*B + D,1e-9);

%% SM FORM
% transfer function must have the same denominator.
den = tf(G.den{1,1},1);
% P = tf(Gp.num,1);
P = zpk(minreal(G*den));

% addpath('/Applications/MATLAB_R2016b.app/toolbox/Mult&T')

% Obtain the smith McMillan form
Psym   = tf2sym(P);
densym = tf2sym(den);
[U_Lsym,U_Rsym,M_Psym] = smithForm(Psym);
M_Psym = M_Psym/densym;
U_Linv = zpk(sym2tf(inv(U_Lsym)));
U_Rinv = zpk(sym2tf(inv(U_Rsym)));
U_L    = zpk(sym2tf(U_Lsym));
U_R    = zpk(sym2tf(U_Rsym));
M_P    = zpk(sym2tf(M_Psym));
M_P    = minreal(M_P);

%% Controller
zeta = 1/sqrt(2);
w0 = 2*pi*5 / sqrt(sqrt(4*zeta^4+1) - 2*zeta^2);

g = w0^2/(s^2+2*zeta*w0*s+w0^2);

M_Y = [1/M_P(1,1) 0; 0 1/M_P(2,2)]*g;
M_T = minreal(M_P*M_Y,1e-4);

Y = minreal(U_R*M_Y*U_L,1e-5);

% Check youla stability
isstable(Y)

T_y = minreal(G*Y,1e-5);
S_y = minreal(eye(2)-T_y,1e-5);
Gc = minreal(Y*(S_y)^-1,1e-4);

w = logspace(-2,3,200);
figure
bodemag(M_T,w,'b')
hold on;
bodemag(eye(2)-M_T,w,'-.r')
bodemag(Y,w,'k')
legend({'T_y','S_y','Y'},'location','best')
saveas(gcf,'Youal_controller_bode.png')

save('Youla_controller.mat')

