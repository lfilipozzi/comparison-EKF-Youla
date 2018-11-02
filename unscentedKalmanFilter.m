classdef unscentedKalmanFilter < matlab.System & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime & ...
        matlab.system.mixin.CustomIcon 
    % unscentedKalmanFilter Simulink implementation of an EKF
    %
    % This files includes a template to implement an unscented Kalman
    % filter in Simulink.
    % The only properties and methods that needs to be modified are:
    % - x_init
    % - nb_state
    % - nb_input
    % - nb_output
    % - stateEquation
    % - outputEquation

    % Public, tunable properties
    properties
        % Number of state of the model
        nb_state  = 5;
        % Number of input of the model
        nb_input  = 1;
        % Number of output of the model
        nb_output = 3;
        % Initial state
        x_init = [0;
                (50/3.6)/0.3250*9.9649;
                (50/3.6)/0.3250;
                (50/3.6)/0.3250;
                (50/3.6)];
        % Sampling time
        Ts = 0.1;
        % Halfshaft compliance
        K_hsf = 2.2918e+03;
        % Halfshaft damping
        b_hsf = 22.9183;
        % Gearbox ratio
        G = 9.9649;
        % Motor inertia
        Jm_R = 0.0351;
        % Wheel radius
        rw = 0.3250;
        % Wheel inertia
        Jw = 0.9;
        % Vehicle mass
        m = 1530;
        % c1
        c1 = 1.2801;    % Burckhardt parameters
        % c2
        c2 = 23.99;
        % c3
        c3 = 0.25;
        % fRLz0
        fRLz0 = 2.9695e+03;
        % fRRz0
        fRRz0 = 2.9695e+03;
        % Cx
        CRx = 4.3302e+04;
        % muRL0
        muRL0 = 1;
        % muRR0
        muRR0 = 1;
        % epsDugoff
        epsDugoff = 0.01;
    end

    properties(DiscreteState)
        x_apriori
        x_aposteriori
        P_apriori
        P_aposteriori
        u_k
        Q_k
    end

    % Pre-computed constants
    properties(Access = private)
        
    end

    methods(Access = protected)
        function setupImpl(~)
            % Perform one-time calculations, such as computing constants
            
        end

        function [x_aposteriori, P_aposteriori] = stepImpl(obj,u_k,y_k,...
                Q_k,R_k)
            % Implement algorithm. Return the a posteriori estimate and the
            % covariance matrix at the following timestep
            
            % Shift register: Save the input and covariance matrix for the
            % next timestep and obtain the one from the last timestep
            u_prev  = obj.u_k;  % u_{k-1}
            obj.u_k = u_k;      % u_{k}
            Q_prev  = obj.Q_k;  % Q_{k-1}
            obj.Q_k = Q_k;      % Q_{k}
            
            % Model update
            % Create sigma-points
            X_sigma_points = ...
                obj.createSigmaPoints(obj.x_aposteriori,obj.P_aposteriori);
            % Propagate each sigma-point through the state equation
            X_sigma_points_prop = zeros(obj.nb_state,2*obj.nb_state);
            for i = 1:2*obj.nb_state
                X_sigma_points_prop(:,i) = ...
                    obj.forwardEuler_stateEquation(X_sigma_points(:,i),...
                    u_prev);
            end
            % Compute a priori estimate
            obj.x_apriori = mean(X_sigma_points_prop,2);
            % obj.P_apriori = cov(X_sigma_points_prop) + Q_prev;
            obj.P_apriori = zeros(obj.nb_state);
            for i = 1:2*obj.nb_state
                obj.P_apriori = obj.P_apriori + ...
                    (X_sigma_points_prop(:,i) - obj.x_apriori) * ...
                    (X_sigma_points_prop(:,i) - obj.x_apriori)';
            end
            obj.P_apriori = obj.P_apriori / (2*obj.nb_state) + Q_prev;
            
            % Measurement update
            % Create sigma-points
            X_sigma_points = ...
                obj.createSigmaPoints(obj.x_apriori,obj.P_apriori);
            % Propagate sigma-point through the output equation
            Y_sigma_points = zeros(obj.nb_output,2*obj.nb_state);
            for i = 1:2*obj.nb_state
                Y_sigma_points(:,i) = ...
                    obj.outputEquation(X_sigma_points(:,i),u_k);
            end
            % Compute mean and variance
            y_mean = mean(Y_sigma_points,2);
            P_XY   = zeros(obj.nb_state,obj.nb_output);
            for i = 1:2*obj.nb_state
                P_XY = P_XY + (X_sigma_points(:,i) - obj.x_apriori) * ...
                    (Y_sigma_points(:,i) - y_mean)';
            end
            P_XY = P_XY / (2*obj.nb_state);
            % P_Y    = cov(Y_sigma_points) + R_k;
            P_Y = zeros(obj.nb_output);
            for i = 1:2*obj.nb_state
                P_Y = P_Y + ...
                    (Y_sigma_points(:,i) - y_mean) * ...
                    (Y_sigma_points(:,i) - y_mean)';
            end
            P_Y = P_Y / (2*obj.nb_state) + R_k;
            % Compute a posteriori estimate
            obj.x_aposteriori = obj.x_apriori + P_XY / P_Y * (y_k-y_mean);
            obj.P_aposteriori = obj.P_apriori - P_XY / P_Y * P_XY';
            
            % Return outputs
            x_aposteriori = obj.x_aposteriori;
            P_aposteriori = obj.P_aposteriori;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.x_aposteriori = obj.x_init;
            obj.P_aposteriori = eye(obj.nb_state)/1000;
            obj.u_k = zeros(obj.nb_input,1);
            obj.Q_k = eye(obj.nb_state);
        end
        
        function sts = getSampleTimeImpl(obj)
            % Define sampling time
            sts = createSampleTime(obj,'Type','Discrete Periodic',...
              'SampleTime',obj.Ts,'OffsetTime',0);
        end
        
        function ds = getDiscreteStateImpl(obj)
            % Return structure of properties with DiscreteState attribute
            ds = obj.x_aposteriori;
        end

        function flag = isInputSizeLockedImpl(~,~)
            % Return true if input size is not allowed to change while
            % system is running
            flag = true;
        end

        function icon = getIconImpl(~)
            % Define icon for System block
            icon = 'Unscented Kalman Filter';
            % icon = {'My','System'}; % Example: multi-line text icon
        end
        
        function [c1, c2] = isOutputFixedSizeImpl(~)
            % Return true if output size is not allowed to change while
            % system is running
            c1 = true;
            c2 = true;
       end

        function [sz_1, sz_2] = getOutputSizeImpl(obj)
            % Return size for each output port
              sz_1 = [obj.nb_state 1]; 
              sz_2 = [obj.nb_state obj.nb_state]; 

            % Example: inherit size from first input port
            % out = propagatedInputSize(obj,1);
        end
        
        function [out1, out2] = getOutputDataTypeImpl(~)
            out1 = 'double';
            out2 = 'double';
        end
        
        function [c1, c2] = isOutputComplexImpl(~)
            c1 = false;
            c2 = false;
        end
        
        function X_sigma_points = createSigmaPoints(obj,mean,covariance)
            % TODO: change code to allow to use Cholesky decomposition of
            % semi-positive definte matrix
            cholesky_dec = chol(obj.nb_state*covariance)';
            
            % Compute sigma points
            X_sigma_points = zeros(obj.nb_state,2*obj.nb_state);
            for i = 1:obj.nb_state
                X_sigma_points(:,i) = mean + cholesky_dec(:,i);
                X_sigma_points(:,obj.nb_state+i) = mean - cholesky_dec(:,i);
            end
        end
        
        function x_next = forwardEuler_stateEquation(obj,x,u_prev)
            % Implement forward-Euler approximation of the state equation
            
            h = obj.Ts;
            x_dot = obj.stateEquation(x,u_prev);
            x_next = x + x_dot * h;
        end
        
        function x_dot = stateEquation(obj,x,u)
            % State equation: returns x_{k+1} as a function of x_k and u_k
            
            % This method is called by the Runge-Kutta method to obtain the
            % state at the next timestep
            
            % Load parameters
            K_hsf = obj.K_hsf;
            b_hsf = obj.b_hsf;
            G     = obj.G;
            Jm_R  = obj.Jm_R;
            rw    = obj.rw;
            Jw    = obj.Jw;
            m     = obj.m;
            c1    = obj.c1;
            c2    = obj.c2;
            c3    = obj.c3;
            fRLz0 = obj.fRLz0;
            fRRz0 = obj.fRRz0;
            CRx   = obj.CRx;
            muRL0 = obj.muRL0;
            muRR0 = obj.muRR0;
            epsDugoff = obj.epsDugoff;
            
            theta_hsf = x(1);
            wm_R      = x(2);
            ww_RL     = x(3);
            ww_RR     = x(4);
            U         = x(5);
            
            tau_m_R = u;
            
            % Tire slip
            sRLx = (rw*ww_RL - U) / U;
            sRRx = (rw*ww_RR - U) / U;

            % % Normal forces (TODO: solve algebraic loop)
            fRLz = fRLz0;% + DFx*ax;
            fRRz = fRRz0;% + DFx*ax;

%             % Tire forces (Burckhardt model)
%             muRLx = c1*(1-exp(-c2*sRLx)) - c3*sRLx;
%             muRRx = c1*(1-exp(-c2*sRRx)) - c3*sRRx;
%             fRLx = muRLx * fRLz;
%             fRRx = muRRx * fRRz;
            
            % Dugoff tire model
            alphaRL = 0;
            alphaRR = 0;
            vRLx = U;
            vRRx = U;
            CRy = 0;
            [fRLx,~] = Dugoff(CRx,CRy,fRLz,sRLx,alphaRL,vRLx,muRL0,epsDugoff);
            [fRRx,~] = Dugoff(CRx,CRy,fRRz,sRRx,alphaRR,vRRx,muRR0,epsDugoff);

            % Halshaft torque
            tau_hsf = K_hsf * theta_hsf + b_hsf * (2/G*wm_R-ww_RL-ww_RR); 

            % Equation of motion
            wm_R_dot      = 1/Jm_R * (tau_m_R - 2/G * tau_hsf);  
            theta_hsf_dot = 2/G * wm_R - ww_RL - ww_RR;
            ww_RL_dot = (tau_hsf - rw * fRLx) /Jw;
            ww_RR_dot = (tau_hsf - rw * fRRx) /Jw;
            U_dot     = 1/m * (fRLx + fRRx);
            
            % Return state derivatives
            x_dot = [theta_hsf_dot wm_R_dot ww_RL_dot ww_RR_dot U_dot]';
        end
        
        function y = outputEquation(obj,x,~)
            % Output equation: returns y_k as a function of x_k and u_k
            
            % Load parameters
            K_hsf = obj.K_hsf;
            b_hsf = obj.b_hsf;
            G     = obj.G;
            Jm_R  = obj.Jm_R;
            rw    = obj.rw;
            Jw    = obj.Jw;
            m     = obj.m;
            c1    = obj.c1;
            c2    = obj.c2;
            c3    = obj.c3;
            fRLz0 = obj.fRLz0;
            fRRz0 = obj.fRRz0;
            CRx   = obj.CRx;
            muRL0 = obj.muRL0;
            muRR0 = obj.muRR0;
            epsDugoff = obj.epsDugoff;
            
            % Compute wRL and wRR
            wRL = x(3);
            wRR = x(4);
            
            % Compute U_dot
            U = x(5);
            % RL longitudinal slip
            vRLx = U;
            if abs(rw*wRL) > abs(vRLx)
                den = rw*wRL;
            else
                den = vRLx;
            end
            if den == 0
                sRLx = 0;
            else
                sRLx = -(vRLx - rw * wRL) / den;
            end
            % RR longitudinal slip
            vRRx = U;
            if abs(rw*wRR) > abs(vRRx)
                den = rw*wRR;
            else
                den = vRRx;
            end
            if den == 0
                sRRx = 0;
            else
                sRRx = -(vRRx - rw * wRR) / den;
            end
%             sRLx  = (rw*wRL - U) / U;
%             sRRx  = (rw*wRR - U) / U;
            
            fRLz  = fRLz0;
            fRRz  = fRRz0;
            
            % Burckhartd tire model
            muRLx = c1*(1-exp(-c2*sRLx)) - c3*sRLx;
            muRRx = c1*(1-exp(-c2*sRRx)) - c3*sRRx;
            fRLx  = muRLx * fRLz;
            fRRx  = muRRx * fRRz;
            
            % Dugoff tire model
            alphaRL = 0;
            alphaRR = 0;
            vRLx = U;
            vRRx = U;
            CRy = 0;
            [fRLx,~] = Dugoff(CRx,CRy,fRLz,sRLx,alphaRL,vRLx,muRL0,epsDugoff);
            [fRRx,~] = Dugoff(CRx,CRy,fRRz,sRRx,alphaRR,vRRx,muRR0,epsDugoff);
            
            % Compute U_dot
            U_dot = 1/m * (fRLx + fRRx);
            
            % Return the output
            y = [U_dot wRL wRR]';
        end
    end
end
