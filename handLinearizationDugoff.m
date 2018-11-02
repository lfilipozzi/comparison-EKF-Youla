function [A,C] = handLinearizationDugoff(ww_RL,ww_RR,U,...
    K_hsf,b_hsf,G,Jm_R,rw,Jw,m,fRLz0,fRRz0,CRx,muRL0,muRR0,epsDugoff)
% x = [theta_hsf wm_R ww_RL ww_RR U];
% u = tau_m_R;
% y = [U_dot ww_RL ww_RR];sRLx = (rw*ww_RL - U) / U;

%% Compute tire forces
% Compute slip ratio
sRLx = (rw*ww_RL - U) / U;
DsRLxDww_RL = rw/U;
DsRLxDww_RR = 0;
DsRLxDU     = -rw*ww_RL/U^2;

sRRx = (rw*ww_RR - U) / U;
DsRRxDww_RL = 0;
DsRRxDww_RR = rw/U;
DsRRxDU     = -rw*ww_RR/U^2;

% RL Dugoff tire model
muRL  = muRL0 * (1 - epsDugoff*U*abs(sRLx));
DmuRLDww_RL = -muRL0 * epsDugoff * U * sign(sRLx) * DsRLxDww_RL;
DmuRLDww_RR = -muRL0 * epsDugoff * U * sign(sRLx) * DsRLxDww_RR;
DmuRLDU     = -muRL0 * epsDugoff * (abs(sRLx) + U * sign(sRLx) * DsRLxDU);

if abs(sRLx) == 1
    % fRLx = sign(sRLx)*(CRx * muRL * fRLz0)/CRx;
    DfRLxDtheta_hsf = 0;
    DfRLxDwm_R      = 0;
    DfRLxDww_RL     = abs(sRLx)*DsRLxDww_RL*(CRx * muRL * fRLz0)/CRx + sign(sRLx)*(CRx * DmuRLDww_RL * fRLz0)/CRx;
    DfRLxDww_RR     = abs(sRLx)*DsRLxDww_RR*(CRx * muRL * fRLz0)/CRx + sign(sRLx)*(CRx * DmuRLDww_RR * fRLz0)/CRx;
    DfRLxDU         = abs(sRLx)*DsRLxDU    *(CRx * muRL * fRLz0)/CRx + sign(sRLx)*(CRx * DmuRLDU     * fRLz0)/CRx;
else
    fRLxd = CRx * sRLx / abs(1-abs(sRLx));
    DfRLxdDww_RL = CRx * (DsRLxDww_RL * abs(1-abs(sRLx)) - sRLx * sign(sRLx*(abs(sRLx)-1))*DsRLxDww_RL) / (1-abs(sRLx))^2;
    DfRLxdDww_RR = CRx * (DsRLxDww_RR * abs(1-abs(sRLx)) - sRLx * sign(sRLx*(abs(sRLx)-1))*DsRLxDww_RR) / (1-abs(sRLx))^2;
    DfRLxdDU     = CRx * (DsRLxDU     * abs(1-abs(sRLx)) - sRLx * sign(sRLx*(abs(sRLx)-1))*DsRLxDU    ) / (1-abs(sRLx))^2;
    
    lambdaRL = muRL / (2 * abs(fRLxd/fRLz0));
    DlambdaRLDww_RL = (DmuRLDww_RL * (2*abs(fRLxd/fRLz0)) - muRL * (2*abs(fRLz0)*sign(fRLxd)*DfRLxdDww_RL)) / (2 * abs(fRLxd/fRLz0))^2;
    DlambdaRLDww_RR = (DmuRLDww_RR * (2*abs(fRLxd/fRLz0)) - muRL * (2*abs(fRLz0)*sign(fRLxd)*DfRLxdDww_RR)) / (2 * abs(fRLxd/fRLz0))^2;
    DlambdaRLDU     = (DmuRLDU     * (2*abs(fRLxd/fRLz0)) - muRL * (2*abs(fRLz0)*sign(fRLxd)*DfRLxdDU    )) / (2 * abs(fRLxd/fRLz0))^2;
    
    if lambdaRL > 1
        % No saturation
        % fRLx = fRLxd;
        DfRLxDtheta_hsf = 0;
        DfRLxDwm_R      = 0;
        DfRLxDww_RL     = DfRLxdDww_RL;
        DfRLxDww_RR     = DfRLxdDww_RR;
        DfRLxDU         = DfRLxdDU;
        
    else
        % Saturation
        % fRLx = fRLxd * (2 * lambdaRL - lambdaRL^2);
        DfRLxDtheta_hsf = 0;
        DfRLxDwm_R      = 0;
        DfRLxDww_RL     = DfRLxdDww_RL * (2 * lambdaRL - lambdaRL^2) + fRLxd * (2 - 2*lambdaRL) * DlambdaRLDww_RL;
        DfRLxDww_RR     = DfRLxdDww_RR * (2 * lambdaRL - lambdaRL^2) + fRLxd * (2 - 2*lambdaRL) * DlambdaRLDww_RR;
        DfRLxDU         = DfRLxdDU     * (2 * lambdaRL - lambdaRL^2) + fRLxd * (2 - 2*lambdaRL) * DlambdaRLDU;
        
    end
end

% RR Dugoff tire model
muRR  = muRR0 * (1 - epsDugoff*U*abs(sRRx));
DmuRRDww_RL = -muRR0 * epsDugoff * U * sign(sRRx) * DsRRxDww_RR;
DmuRRDww_RR = -muRR0 * epsDugoff * U * sign(sRRx) * DsRRxDww_RR;
DmuRRDU     = -muRR0 * epsDugoff * (abs(sRRx) + U * sign(sRRx) * DsRRxDU);

if abs(sRRx) == 1
    % fRRx = sign(sRRx)*(CRx * muRR * fRRz0)/CRx;
    DfRRxDtheta_hsf = 0;
    DfRRxDwm_R      = 0;
    DfRRxDww_RL     = abs(sRRx)*DsRRxDww_RL*(CRx * muRR * fRRz0)/CRx + sign(sRRx)*(CRx * DmuRRDww_RL * fRRz0)/CRx;
    DfRRxDww_RR     = abs(sRRx)*DsRRxDww_RR*(CRx * muRR * fRRz0)/CRx + sign(sRRx)*(CRx * DmuRRDww_RR * fRRz0)/CRx;
    DfRRxDU         = abs(sRRx)*DsRRxDU    *(CRx * muRR * fRRz0)/CRx + sign(sRRx)*(CRx * DmuRRDU     * fRRz0)/CRx;
else
    fRRxd = CRx * sRRx / abs(1-abs(sRRx));
    DfRRxdDww_RL = CRx * (DsRRxDww_RL * abs(1-abs(sRRx)) - sRRx * sign(sRRx*(abs(sRRx)-1))*DsRRxDww_RL) / (1-abs(sRRx))^2;
    DfRRxdDww_RR = CRx * (DsRRxDww_RR * abs(1-abs(sRRx)) - sRRx * sign(sRRx*(abs(sRRx)-1))*DsRRxDww_RR) / (1-abs(sRRx))^2;
    DfRRxdDU     = CRx * (DsRRxDU     * abs(1-abs(sRRx)) - sRRx * sign(sRRx*(abs(sRRx)-1))*DsRRxDU    ) / (1-abs(sRRx))^2;
    
    lambdaRR = muRR / (2 * abs(fRRxd/fRRz0));
    DlambdaRRDww_RL = (DmuRRDww_RL * (2*abs(fRRxd/fRRz0)) - muRR * (2*abs(fRRz0)*sign(fRRxd)*DfRRxdDww_RL)) / (2 * abs(fRRxd/fRRz0))^2;
    DlambdaRRDww_RR = (DmuRRDww_RR * (2*abs(fRRxd/fRRz0)) - muRR * (2*abs(fRRz0)*sign(fRRxd)*DfRRxdDww_RR)) / (2 * abs(fRRxd/fRRz0))^2;
    DlambdaRRDU     = (DmuRRDU     * (2*abs(fRRxd/fRRz0)) - muRR * (2*abs(fRRz0)*sign(fRRxd)*DfRRxdDU    )) / (2 * abs(fRRxd/fRRz0))^2;
    
    if lambdaRR > 1
        % No saturation
        % fRRx = fRRxd;
        DfRRxDtheta_hsf = 0;
        DfRRxDwm_R      = 0;
        DfRRxDww_RL     = DfRRxdDww_RL;
        DfRRxDww_RR     = DfRRxdDww_RR;
        DfRRxDU         = DfRRxdDU;
        
    else
        % Saturation
        % fRRx = fRRxd * (2 * lambdaRR - lambdaRR^2);
        DfRRxDtheta_hsf = 0;
        DfRRxDwm_R      = 0;
        DfRRxDww_RL     = DfRRxdDww_RL * (2 * lambdaRR - lambdaRR^2) + fRRxd * (2 - 2*lambdaRR) * DlambdaRRDww_RL;
        DfRRxDww_RR     = DfRRxdDww_RR * (2 * lambdaRR - lambdaRR^2) + fRRxd * (2 - 2*lambdaRR) * DlambdaRRDww_RR;
        DfRRxDU         = DfRRxdDU     * (2 * lambdaRR - lambdaRR^2) + fRRxd * (2 - 2*lambdaRR) * DlambdaRRDU;
        
    end
end

%% Compute halfshaft torque
% tau_hsf = K_hsf * theta_hsf + b_hsf * (2/G*wm_R - ww_RL - ww_RR); 
Dtau_hsfDtheta_hsf = K_hsf;
Dtau_hsfDwm_R      = b_hsf*2/G;
Dtau_hsfDww_RL     = -b_hsf;
Dtau_hsfDww_RR     = -b_hsf;
Dtau_hsfDU         = 0;

%% Equation of motion
% wm_R_dot      = 1/Jm_R * (tau_m_R - 2/G * tau_hsf);  
Dwm_R_dotDtheta_hsf = 1/Jm_R * (- 2/G * Dtau_hsfDtheta_hsf); 
Dwm_R_dotDwm_R      = 1/Jm_R * (- 2/G * Dtau_hsfDwm_R     ); 
Dwm_R_dotDww_RL     = 1/Jm_R * (- 2/G * Dtau_hsfDww_RL    ); 
Dwm_R_dotDww_RR     = 1/Jm_R * (- 2/G * Dtau_hsfDww_RR    ); 
Dwm_R_dotDU         = 1/Jm_R * (- 2/G * Dtau_hsfDU        ); 

% theta_hsf_dot = 2/G * wm_R - ww_RL - ww_RR;
Dtheta_hsf_dotDtheta_hsf = 0;
Dtheta_hsf_dotDwm_R      = 2/G;
Dtheta_hsf_dotDww_RL     = -1;
Dtheta_hsf_dotDww_RR     = -1;
Dtheta_hsf_dotDU         = 0;

% ww_RL_dot = (tau_hsf - rw * fRLx) /Jw;
Dww_RL_dotDtheta_hsf = (Dtau_hsfDtheta_hsf - rw * DfRLxDtheta_hsf) /Jw;
Dww_RL_dotDwm_R      = (Dtau_hsfDwm_R      - rw * DfRLxDwm_R     ) /Jw;
Dww_RL_dotDww_RL     = (Dtau_hsfDww_RL     - rw * DfRLxDww_RL    ) /Jw;
Dww_RL_dotDww_RR     = (Dtau_hsfDww_RR     - rw * DfRLxDww_RR    ) /Jw;
Dww_RL_dotDU         = (Dtau_hsfDU         - rw * DfRLxDU        ) /Jw;

% ww_RR_dot = (tau_hsf - rw * fRRx) /Jw;
Dww_RR_dotDtheta_hsf = (Dtau_hsfDtheta_hsf - rw * DfRRxDtheta_hsf) /Jw;
Dww_RR_dotDwm_R      = (Dtau_hsfDwm_R      - rw * DfRRxDwm_R     ) /Jw;
Dww_RR_dotDww_RL     = (Dtau_hsfDww_RL     - rw * DfRRxDww_RL    ) /Jw;
Dww_RR_dotDww_RR     = (Dtau_hsfDww_RR     - rw * DfRRxDww_RR    ) /Jw;
Dww_RR_dotDU         = (Dtau_hsfDU         - rw * DfRRxDU        ) /Jw;

% U_dot     = 1/m * (fRLx + fRRx);
DU_dotDtheta_hsf = 1/m * (DfRLxDtheta_hsf + DfRRxDtheta_hsf);
DU_dotDwm_R      = 1/m * (DfRLxDwm_R      + DfRRxDwm_R     );
DU_dotDww_RL     = 1/m * (DfRLxDww_RL     + DfRRxDww_RL    );
DU_dotDww_RR     = 1/m * (DfRLxDww_RR     + DfRRxDww_RR    );
DU_dotDU         = 1/m * (DfRLxDU         + DfRRxDU        );

%% Compute Jacobian matrices A and C
% x = [theta_hsf wm_R ww_RL ww_RR U];
A = [Dtheta_hsf_dotDtheta_hsf Dtheta_hsf_dotDwm_R Dtheta_hsf_dotDww_RL Dtheta_hsf_dotDww_RR Dtheta_hsf_dotDU;
     Dwm_R_dotDtheta_hsf      Dwm_R_dotDwm_R      Dwm_R_dotDww_RL      Dwm_R_dotDww_RR      Dwm_R_dotDU;
     Dww_RL_dotDtheta_hsf     Dww_RL_dotDwm_R     Dww_RL_dotDww_RL     Dww_RL_dotDww_RR     Dww_RL_dotDU;
     Dww_RR_dotDtheta_hsf     Dww_RR_dotDwm_R     Dww_RR_dotDww_RL     Dww_RR_dotDww_RR     Dww_RR_dotDU;
     DU_dotDtheta_hsf         DU_dotDwm_R         DU_dotDww_RL         DU_dotDww_RR         DU_dotDU];

% y = [U_dot ww_RL ww_RR];
C = [DU_dotDtheta_hsf DU_dotDwm_R DU_dotDww_RL DU_dotDww_RR DU_dotDU;
     0                0           1            0            0;
     0                0           0            1            0];

end







