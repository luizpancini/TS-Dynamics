function [tp,d,x0_save,inputs_outputs] = AEBL_rkf_Sheng_Lyapunov(tspan,x0,xdot0,inputs_outputs,d_i,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options)

%% Handle the inputs
if ~isempty(options)
    if any(strcmp(options,'hlim')) == 1, k=find(strcmp(options,'hlim')); hlim = options{k+1}; end
    if any(strcmp(options,'xdotRelTol')) == 1, k=find(strcmp(options,'xdotRelTol')); xdotRelTol = options{k+1}; end
    if any(strcmp(options,'xdotAbsTol')) == 1, k=find(strcmp(options,'xdotAbsTol')); xdotAbsTol = options{k+1}; end
    if any(strcmp(options,'relax_fac_min')) == 1, k=find(strcmp(options,'relax_fac_min')); relax_fac_min = options{k+1}; end
    if any(strcmp(options,'alpha_bar_lim')) == 1, k=find(strcmp(options,'alpha_bar_lim')); alpha_bar_lim = options{k+1}; end
    if any(strcmp(options,'xdot_max_it')) == 1, k=find(strcmp(options,'xdot_max_it')); xdot_max_it = options{k+1}; end
    if any(strcmp(options,'FD_method')) == 1, k=find(strcmp(options,'FD_method')); FD_method = options{k+1}; end
    if any(strcmp(options,'int_method')) == 1, k=find(strcmp(options,'int_method')); int_method = options{k+1}; end
    if any(strcmp(options,'d0')) == 1, k=find(strcmp(options,'d0')); d0 = options{k+1}; end
    if any(strcmp(options,'r')) == 1, k=find(strcmp(options,'r')); r = options{k+1}; end
end        

if ~exist('hlim','var'), hlim = 1e-5; end                               % Default maximum timestep
if ~exist('xdotRelTol','var'), xdotRelTol = 1e-12; end                  % Default xdot relative tolerance
if ~exist('xdotAbsTol','var'), xdotAbsTol = 1e-12; end                  % Default xdot absolute tolerance
if ~exist('relax_fac_min','var'), relax_fac_min = 2^-3; end             % Default minimum relaxation factor (0-1): decrease for better (but slower) convergence on xdot
if ~exist('alpha_bar_lim','var'), alpha_bar_lim = pi/4; end             % Default maximum effective AoA for DS model validity
if ~exist('xdot_max_it','var'), xdot_max_it = 10; end                   % Default maximum number of iterations for xdot implicit convergence
if ~exist('FD_method','var'), FD_method = 'c2'; end                     % Default finite-difference method (2-point centered)
if ~exist('int_method','var'), int_method = 1; end                      % Default integration method of Jacobian (1 for Euler, 2 for Ralston, 3 for RK4)
if ~exist('d0','var'), d0 = 1; end                                      % Default disturbance norm
if ~exist('r','var'), r = 1; end                                        % Default time steps for renormalization

%% Setup the algorithm
% Time variables
ti = tspan(1);          % Initial time
tf = tspan(end);        % Final time
tc = ti;                % Current time
% Pre-allocate
s0 = round((tf-ti)/hlim)+1;  % Initial size of the vectors
tp = zeros(s0,1); d = zeros(length(x0),s0);
% Set initial conditions
x_i = x0;
xdot_g = xdot0;
% Initialize outputs
d(:,1) = d_i;           
% Initialize global variables
global divergent_flag; divergent_flag = 0;
% Jacobian calculation variables
N = length(x0);       % Jacobian size
dx = 1e-10;           % Finite difference

%% Solve the ODEs
i = 1;
hc = hlim;                          % Set fixed timestep 
while tc < tf
    relax_fac = 1;                  % Reset relaxation factor
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
    % Check for divergent oscillations
    if abs(inputs_outputs(4,2)) > alpha_bar_lim
        warning(['Effective AoA > ', num2str(alpha_bar_lim*180/pi),'° - stopping integration'])
        divergent_flag = 1;
        break
    end
    % Calculate linearized disturbance
    J = AEBL_Sheng_Jacobian_int(int_method,FD_method,N,tc,x_i,xdot_g,inputs_outputs,k1,x2,x3,x4,x5,hc,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,dx);
    wdot = J*d_i;
    d_ip1 = d_i + wdot*hc;
    % Save outputs   
    d(:,i+1) = d_ip1;
    % Renormalization
    if rem(i+1,r) == 0
        d_ip1 = d0*d_ip1/norm(d_ip1); 
    end
    % Setup next time step 
    tc = tc+hc;
    tp(i+1) = tc;
    x_i = x_ip1;
    d_i = d_ip1;
    i = i+1;
    % Display progress
    if rem(i,1e4) == 0
        disp(['RKF45 progress: ',num2str((tc-ti)/(tf-ti)*100,'%10.2f') '%'])
    end
end
x0_save = x_i;
            
% Truncate pre-allocated vectors
tp = tp(1:s0); d = d(:,1:s0);
