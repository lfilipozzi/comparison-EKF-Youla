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