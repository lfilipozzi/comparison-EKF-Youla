classdef extendedKalmanFilter < matlab.System & ...
        matlab.system.mixin.Propagates & ...
        matlab.system.mixin.SampleTime & ...
        matlab.system.mixin.CustomIcon 
    % extendedKalmanFilter Simulink implementation of an EKF
    %
    % This files includes a template to implement an extended Kalman filter
    % in Simulink.
    % The only properties and methods that needs to be modified are:
    % - x_init
    % - nb_state
    % - nb_input
    % - nb_output
    % - linearizedStateMatrices
    % - linearizedOutputMatrices
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
        % Burckhardt parameters
        c1 = 1.2801;
        c2 = 23.99;
        c3 = 0.25;
        % Normal forces
        fRLz0 = 2.9695e+03;
        fRRz0 = 2.9695e+03;
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
        function setupImpl(obj)
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
            [A,E] = obj.linearizedStateMatrices();
            %obj.x_apriori = stateEquation(obj.x_aposteriori,u_prev);
            %obj.x_apriori = rungeKutta4_stateEquation(obj,u_prev,u_k,h);
            obj.x_apriori = obj.forwardEuler_stateEquation(u_prev);
            obj.P_apriori = A*obj.P_aposteriori*A' + E*Q_prev*E';
            
            % Measurement update
            [C,F] = obj.linearizedOutputMatrices();
            L = obj.P_apriori*C' / (C*obj.P_apriori*C' + F*R_k*F');
            obj.x_aposteriori = obj.x_apriori + L * (y_k - ...
                obj.outputEquation(obj.x_apriori,u_k));
            obj.P_aposteriori = obj.P_apriori - L*C*obj.P_apriori;
            
            % Return outputs
            x_aposteriori = obj.x_aposteriori;
            P_aposteriori = obj.P_aposteriori;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.x_aposteriori = obj.x_init;
            obj.P_aposteriori = zeros(obj.nb_state);
            obj.u_k = zeros(obj.nb_input,1);
            obj.Q_k = zeros(obj.nb_input);
        end
        
        function sts = getSampleTimeImpl(obj)
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
            icon = 'Extended Kalman Filter';
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
        
%         function x_next = rungeKutta4_stateEquation(obj,u_prev,u_k)
%             % Implement 4-th order Runge-Kutta method of the state equation
%             
%             % TODO: Why using RK4 when Forward-Euler is used to compute the A and
%             % C matrices? Replace to use Forward Euler
%             
%             x = obj.x_aposteriori;
%             h = obj.Ts;
%             
%             % TODO: obtain the value of u2 (the value of u1 has already
%             % been saved in the discrete state)
%             u1 = u_prev;    % input at t = tn
%             u2 = % input at t = tn + h/2
%             u3 = u_k;       % input at t = tn + h
%             
%             k1 = h * obj.stateEquation(x     ,u1);
%             k2 = h * obj.stateEquation(x+k1/2,u2);
%             k3 = h * obj.stateEquation(x+k2/2,u2);
%             k4 = h * obj.stateEquation(x+k3  ,u3);
%             
%             x_next = x + (k1 + 2*k2 + 2*k3 + k4)/6;
%         end
        
        function x_next = forwardEuler_stateEquation(obj,u_prev)
            % Implement forward-Euler approximation of the state equation
            
            x = obj.x_aposteriori;
            h = obj.Ts;
            
            x_dot = obj.stateEquation(x,u_prev);
            
            x_next = x + x_dot * h;
        end
        
        function [A,E] = linearizedStateMatrices(obj)
            % Return the A and E matrices.
            % A is the linearized state equation with respect to the state
            % evaluated at the previous a posteriori estimate.
            % E is the linearized state equation with respect to the noise
            % evaluated at the previous a posteriori estimate.
            
            x = obj.x_aposteriori;
            
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
            
            theta_hsf = x(1);
            wm_R  = x(2);
            ww_RL = x(3);
            ww_RR = x(4);
            U     = x(5);
            
            A = [                 0,                   2/G,                                                                           -1,                                                                           -1,                                                                                                                                                                                                                      0;
                -(2*K_hsf)/(G*Jm_R), -(4*b_hsf)/(G^2*Jm_R),                                                           (2*b_hsf)/(G*Jm_R),                                                           (2*b_hsf)/(G*Jm_R),                                                                                                                                                                                                                      0;
                           K_hsf/Jw,      (2*b_hsf)/(G*Jw), -(b_hsf - fRLz0*rw*((c3*rw)/U - (c1*c2*rw*exp((c2*(U - rw*ww_RL))/U))/U))/Jw,                                                                    -b_hsf/Jw,                                                                                                        (fRLz0*rw*((c3*(U - rw*ww_RL))/U^2 - c3/U + c1*exp((c2*(U - rw*ww_RL))/U)*(c2/U - (c2*(U - rw*ww_RL))/U^2)))/Jw;
                           K_hsf/Jw,      (2*b_hsf)/(G*Jw),                                                                    -b_hsf/Jw, -(b_hsf - fRRz0*rw*((c3*rw)/U - (c1*c2*rw*exp((c2*(U - rw*ww_RR))/U))/U))/Jw,                                                                                                        (fRRz0*rw*((c3*(U - rw*ww_RR))/U^2 - c3/U + c1*exp((c2*(U - rw*ww_RR))/U)*(c2/U - (c2*(U - rw*ww_RR))/U^2)))/Jw;
                                  0,                     0,             -(fRLz0*((c3*rw)/U - (c1*c2*rw*exp((c2*(U - rw*ww_RL))/U))/U))/m,             -(fRRz0*((c3*rw)/U - (c1*c2*rw*exp((c2*(U - rw*ww_RR))/U))/U))/m, -(fRLz0*((c3*(U - rw*ww_RL))/U^2 - c3/U + c1*exp((c2*(U - rw*ww_RL))/U)*(c2/U - (c2*(U - rw*ww_RL))/U^2)) + fRRz0*((c3*(U - rw*ww_RR))/U^2 - c3/U + c1*exp((c2*(U - rw*ww_RR))/U)*(c2/U - (c2*(U - rw*ww_RR))/U^2)))/m];
            A = eye(obj.nb_state) + A * obj.Ts;
            E = eye(obj.nb_input);
        end
        
        function [C,F] = linearizedOutputMatrices(obj)
            % Return the C, and F matrices.
            % C is the linearized output equation with respect to the state
            % evaluated at the current a priori estimate.
            % F is the linearized output equation with respect to the noise
            % evaluated at the current a priori estimate.
            
            x = obj.x_apriori;
            
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
            
            % theta_hsf = x(1);
            % wm_R  = x(2);
            ww_RL = x(3);
            ww_RR = x(4);
            U     = x(5);
            
            C13 = -(fRLz0*((c3*rw)/U-(c1*c2*rw*exp((c2*(U-rw*ww_RL))/U))...
                /U))/m;
            C14 = -(fRRz0*((c3*rw)/U-(c1*c2*rw*exp((c2*(U-rw*ww_RR))/U))...
                /U))/m;
            C15 = -(fRLz0*((c3*(U-rw*ww_RL))/U^2-c3/U+c1*...
                exp((c2*(U-rw*ww_RL))/U)*(c2/U-(c2*(U-rw*ww_RL))/U^2))+...
                fRRz0*((c3*(U-rw*ww_RR))/U^2-c3/U+c1*...
                exp((c2*(U-rw*ww_RR))/U)*(c2/U-(c2*(U-rw*ww_RR))/U^2)))/m;
            C = [0 0 C13 C14 C15;
                 0 0 1   0   0;
                 0 0 0   1   0];
            F = eye(obj.nb_output);
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

            % Tire forces 5Burckhardt model)
            muRLx = c1*(1-exp(-c2*sRLx)) - c3*sRLx;
            muRRx = c1*(1-exp(-c2*sRRx)) - c3*sRRx;

            fRLx = muRLx * fRLz;
            fRRx = muRRx * fRRz;

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
            
            % Compute wRL and wRR
            wRL = x(3);
            wRR = x(4);
            
            % Compute U_dot
            U = x(5);
            sRLx  = (rw*wRL - U) / U;
            sRRx  = (rw*wRR - U) / U;
            fRLz  = fRLz0;
            fRRz  = fRRz0;
            muRLx = c1*(1-exp(-c2*sRLx)) - c3*sRLx;
            muRRx = c1*(1-exp(-c2*sRRx)) - c3*sRRx;
            fRLx  = muRLx * fRLz;
            fRRx  = muRRx * fRRz;
            U_dot = 1/m * (fRLx + fRRx);
            
            % Return the output
            y = [U_dot wRL wRR]';
        end
    end
end
