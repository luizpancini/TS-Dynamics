function [tp,w,x_i,inputs_outputs] = AEBL_rkf_Sheng_Lyapunov_2PM(tspan,x0,xdot0,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options)

%% Handle the inputs
if ~isempty(options)
    if any(strcmp(options,'hlim')) == 1, k=find(strcmp(options,'hlim')); hlim = options{k+1}; end
    if any(strcmp(options,'xdotRelTol')) == 1, k=find(strcmp(options,'xdotRelTol')); xdotRelTol = options{k+1}; end
    if any(strcmp(options,'xdotAbsTol')) == 1, k=find(strcmp(options,'xdotAbsTol')); xdotAbsTol = options{k+1}; end
    if any(strcmp(options,'relax_fac_min')) == 1, k=find(strcmp(options,'relax_fac_min')); relax_fac_min = options{k+1}; end
    if any(strcmp(options,'alpha_bar_lim')) == 1, k=find(strcmp(options,'alpha_bar_lim')); alpha_bar_lim = options{k+1}; end
    if any(strcmp(options,'xdot_max_it')) == 1, k=find(strcmp(options,'xdot_max_it')); xdot_max_it = options{k+1}; end
    if any(strcmp(options,'d0')) == 1, k=find(strcmp(options,'d0')); d0 = options{k+1}; end
    if any(strcmp(options,'r')) == 1, k=find(strcmp(options,'r')); r = options{k+1}; end
end        

if ~exist('hlim','var'), hlim = 1e-5; end                               % Default maximum timestep
if ~exist('xdotRelTol','var'), xdotRelTol = 1e-12; end                  % Default xdot relative tolerance
if ~exist('xdotAbsTol','var'), xdotAbsTol = 1e-12; end                  % Default xdot absolute tolerance
if ~exist('relax_fac_min','var'), relax_fac_min = 2^-3; end             % Default minimum relaxation factor (0-1): decrease for better (but slower) convergence on xdot
if ~exist('alpha_bar_lim','var'), alpha_bar_lim = pi/4; end             % Default maximum effective AoA for DS model validity
if ~exist('xdot_max_it','var'), xdot_max_it = 10; end                   % Default maximum number of iterations for xdot implicit convergence
if ~exist('d0','var'), d0 = 1e-9; end                                   % Default disturbance norm
if ~exist('r','var'), r = 1; end                                        % Default time steps for renormalization

%% Setup the algorithm
% Time variables
ti = tspan(1);          % Initial time
tf = tspan(end);        % Final time
tc = ti;                % Current time
% Pre-allocate
s0 = round((tf-ti)/hlim+1);  % Initial size of the vectors
tp = zeros(s0,1); w = zeros(s0,length(x0));
% Initialize Lyapunov exponent variables
w0 = ones(1,length(x0)); w0 = w0/norm(w0)*d0;   % Unit disturbance vector
w(1,:) = w0;                                    % Difference vector
% Set initial conditions
x_i = x0; xd_i = x0+w0;
xdot_g = xdot0; xdot_gd = xdot0; 
inputs_outputsd = inputs_outputs;
% Initialize outputs
tp(1) = ti;      
global divergent_flag;

%% Solve the ODEs
i = 1;
hc = hlim;                % Reset current timestep to maximum allowable
k = 1;
while tp(i) < tf
    % Calculate undisturbed flow 
    relax_fac = 1;            % Reset relaxation factor
    if tc+hc > tf, hc = tf-tc; end  % Adjust for last time step
    xdot_it = 0;                    % Reset xdot iteration
    xdot_eps_rel = 10*xdotRelTol;   % Reset relative xdot epsilon
    xdot_eps_abs = 10*xdotAbsTol;   % Reset absolute xdot epsilon
    while xdot_eps_rel > xdotRelTol || xdot_eps_abs > xdotAbsTol
        [k1,inputs_outputs] = AEBL_dxdt_Sheng(tc,x_i,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x1 = x_i + k1*hc/4;
        [k2,~] = AEBL_dxdt_Sheng(tc+hc/4,x1,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x2 = x_i + (3*k1/32+9/32*k2)*hc;
        [k3,~] = AEBL_dxdt_Sheng(tc+3/8*hc,x2,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x3 = x_i + (1932*k1-7200*k2+7296*k3)/2197*hc;
        [k4,~] = AEBL_dxdt_Sheng(tc+12/13*hc,x3,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x4 = x_i + (439/216*k1-8*k2+3680/513*k3-845/4104*k4)*hc;
        [k5,inputs_outputs] = AEBL_dxdt_Sheng(tc+hc,x4,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x5 = x_i + (-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)*hc;
        [k6,~] = AEBL_dxdt_Sheng(tc+hc/2,x5,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        xdot_ng = (16/135*k1+6656/12825*k3+28561/56430*k4-9/50*k5+2/55*k6); 
        xdot_eps_rel = max(abs((xdot_ng-xdot_g)./xdot_g));
        xdot_eps_abs_new = max(abs((xdot_ng-xdot_g)));
        if xdot_eps_abs_new >= xdot_eps_abs 
            if relax_fac>=2*relax_fac_min
                relax_fac = relax_fac/2;
            else
                relax_fac = relax_fac_min;
            end
        end
        xdot_eps_abs = xdot_eps_abs_new;
        xdot_g = relax_fac*xdot_ng+(1-relax_fac)*xdot_g; 
        xdot_it = xdot_it+1;
        if xdot_it==xdot_max_it
            break
        end
    end
    x_ip1 = x_i + xdot_g*hc;
    % Calculate disturbed flow
    relax_fac = 1;            % Reset relaxation factor
    if tc+hc > tf, hc = tf-tc; end  % Adjust for last time step
    xdot_it = 0;                    % Reset xdot iteration
    xdot_eps_rel = 10*xdotRelTol;   % Reset relative xdot epsilon
    xdot_eps_abs = 10*xdotAbsTol;   % Reset absolute xdot epsilon
    while xdot_eps_rel > xdotRelTol || xdot_eps_abs > xdotAbsTol
        [k1d,inputs_outputsd] = AEBL_dxdt_Sheng(tc,xd_i,xdot_gd,tc,inputs_outputsd,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x1d = xd_i + k1d*hc/4;
        [k2d,~] = AEBL_dxdt_Sheng(tc+hc/4,x1d,xdot_gd,tc,inputs_outputsd,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x2d = xd_i + (3*k1d/32+9/32*k2d)*hc;
        [k3d,~] = AEBL_dxdt_Sheng(tc+3/8*hc,x2d,xdot_gd,tc,inputs_outputsd,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x3d = xd_i + (1932*k1d-7200*k2d+7296*k3d)/2197*hc;
        [k4d,~] = AEBL_dxdt_Sheng(tc+12/13*hc,x3d,xdot_gd,tc,inputs_outputsd,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x4d = xd_i + (439/216*k1d-8*k2d+3680/513*k3d-845/4104*k4d)*hc;
        [k5d,inputs_outputsd] = AEBL_dxdt_Sheng(tc+hc,x4d,xdot_gd,tc,inputs_outputsd,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x5d = xd_i + (-8/27*k1d+2*k2d-3544/2565*k3d+1859/4104*k4d-11/40*k5d)*hc;
        [k6d,~] = AEBL_dxdt_Sheng(tc+hc/2,x5d,xdot_gd,tc,inputs_outputsd,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        xdot_ngd = (16/135*k1d+6656/12825*k3d+28561/56430*k4d-9/50*k5d+2/55*k6d); 
        xdot_eps_rel = max(abs((xdot_ngd-xdot_gd)./xdot_gd));
        xdot_eps_abs_new = max(abs((xdot_ngd-xdot_gd)));
        if xdot_eps_abs_new >= xdot_eps_abs 
            if relax_fac>=2*relax_fac_min
                relax_fac = relax_fac/2;
            else
                relax_fac = relax_fac_min;
            end
        end
        xdot_eps_abs = xdot_eps_abs_new;
        xdot_gd = relax_fac*xdot_ngd+(1-relax_fac)*xdot_gd; 
        xdot_it = xdot_it+1;
        if xdot_it == xdot_max_it
            break
        end
    end
    xd_ip1 = xd_i + xdot_gd*hc;
    % Check for divergent oscillations
    if abs(inputs_outputs(4,2)) > alpha_bar_lim
        warning(['Effective AoA > ', num2str(alpha_bar_lim*180/pi),'° - stopping integration'])
        divergent_flag = 1;
        break
    end
    % Calculate difference
    w(i+1,:) = xd_ip1 - x_ip1;
    % Renormalization
    if round(rem(i+1,r)) == 0
        % Renormalize disturbance
        xd_ip1 = x_ip1+w0;
        xdot_gd = xdot_g;
        inputs_outputsd = inputs_outputs;
        % Increase renormalizations count
        k = k+1;
    end
    % Setup next time step
    tc = tc+hc;
    tp(i+1) = tc;
    x_i = x_ip1;
    xd_i = xd_ip1;
    i=i+1;
    if rem(i,1e4) == 0
        disp(['RKF45 progress: ',num2str((tc-ti)/(tf-ti)*100,'%10.2f') '%'])
    end
end

% Truncate pre-allocated vectors
tp = tp(1:s0); w = w(1:s0,:); 
