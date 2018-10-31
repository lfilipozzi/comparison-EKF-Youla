function [mux, muy] = magicFormula(sx,sy,B,C,D)
%MAGICFORMULA Compute friction coef using Pacejka's simplified formula
% This functions requires the longitudinal and lateral slip ratio defined
% as:
%        rw * omega - vx
% sx = -------------------
%      max(vx, rw * omega)
%
%                    vy
% sy = - (1 + sx) * ----
%                   |vx|
% where vx and vy are the longitudinal and lateral velocities of the wheel
% and omega is the wheel angular velocity.

s = sqrt(sx^2 + sy^2);
mu = D*sin(C*atan(B*s));

if s == 0
    mux = 0;
    muy = 0;
else
    mux = sx * mu / s;
    muy = sy * mu / s;
end
    
end

