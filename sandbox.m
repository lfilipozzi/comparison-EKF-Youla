w = .7*2*pi;
tL = 0.5;
A = 25*pi/180;

dt = 0.001/1000;
tmax = 7;
t = 0:dt:tmax;

y = zeros(size(t));

T = 2*pi/w;
i1 = find(t > 3/4*T, 1);
i2 = find(t > 3/4*T + tL, 1);
i3 = find(t > 3/4*T + tL + 1/4*T, 1);
y(1:i1) = A*sin(w*t(1:i1));
y((i1+1):i2) = -A*ones(i2-i1,1);
y((i2+1):i3) = -A*cos(w*(t((i2+1):i3) - t(i2+1)));

plot(t,y*180/pi)

% SWA_swdwell = [t;y]';

% save sinewdwell_SWA SWA_swdwell

%% Plot Dugoff
alpha = 0;
s = linspace(-2,2,100);
fz = m*g/4; 
Cx = CRx;
Cy = CRy;
U = 50/3.6;
mu0 = .5;

fx = zeros(size(s));
for i = 1:length(s)
    [fx(i),~] = Dugoff(Cx,Cy,fz,s(i),alpha,U,mu0,epsDugoff);
end

figure(10)
hold on
plot(s,fx)
    

%% Check observability
U = 80/3.6;
ratio = .9;
ww_RL = ratio*U/rw;
ww_RR = ratio*U/rw;
[A,C] = handLinearizationDugoff(ww_RL,ww_RR,U,...
    K_hsf,b_hsf,G,Jm_R,rw,Jw,m,fRLz0,fRRz0,CRx,muRL0,muRR0,epsDugoff);
rank(obsv(A,C))
    
    