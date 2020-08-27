function [tp,xp,xdot_0,Psi_int] = AEBL_RKF45_Shooting_and_Psi(tspan,x0,xdot0g,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options,dx,Psi_method)

%% Handle the inputs
if ~isempty(options)
    if any(strcmp(options,'hlim')) == 1, k=find(strcmp(options,'hlim')); hlim = options{k+1}; end
    if any(strcmp(options,'RKFtol')) == 1, k=find(strcmp(options,'RKFtol')); RKFtol = options{k+1}; end
    if any(strcmp(options,'xdotRelTol')) == 1, k=find(strcmp(options,'xdotRelTol')); xdotRelTol = options{k+1}; end
    if any(strcmp(options,'xdotAbsTol')) == 1, k=find(strcmp(options,'xdotAbsTol')); xdotAbsTol = options{k+1}; end
    if any(strcmp(options,'relax_fac_min')) == 1, k=find(strcmp(options,'relax_fac_min')); relax_fac_min = options{k+1}; end
    if any(strcmp(options,'xdot_max_it')) == 1, k=find(strcmp(options,'xdot_max_it')); xdot_max_it = options{k+1}; end
    if any(strcmp(options,'b_max_it')) == 1, k=find(strcmp(options,'b_max_it')); b_max_it = options{k+1}; end
    if any(strcmp(options,'RKF_it_max')) == 1, k=find(strcmp(options,'RKF_it_max')); RKF_it_max = options{k+1}; end
end        

if ~exist('hlim','var'), hlim = 1e-4; end                               % Default maximum timestep
if ~exist('RKFtol','var'), RKFtol = 1e-10; end                          % Default RFK tolerance
if ~exist('xdotRelTol','var'), xdotRelTol = 1e-12; end                  % Default xdot relative tolerance
if ~exist('xdotAbsTol','var'), xdotAbsTol = 1e-12; end                  % Default xdot absolute tolerance
if ~exist('relax_fac_min','var'), relax_fac_min = 1; end                % Default minimum relaxation factor (0-1): decrease for better (but slower) convergence on xdot
if ~exist('xdot_max_it','var'), xdot_max_it = 10; end                   % Default maximum number of iterations for xdot implicit convergence
if ~exist('RKF_it_max','var'), RKF_it_max = 10; end                     % Default maximum number of iterations for RKF approximations
if ~exist('b_max_it','var'), b_max_it = 10; end                         % Default maximum number of iterations on boundary

%% Setup the algorithm
% Boundary
delta_b = -1e-16;       % Boundary tolerance
boundary = 1;           % Initialize boundary 
% Time variables
ti = tspan(1);          % Initial time
tf = tspan(end);        % Final time
tc = ti;                % Current time
% Pre-allocate
p_ratio = 1+sqrt(tf-ti);           % Pre-allocation ratio relative to constant time step - my simple formula
s0 = round(p_ratio*(tf-ti)/hlim);  % Initial size of the vectors
tp = zeros(s0,1); xp = zeros(s0,length(x0)); x = xp;
% Set initial conditions
x(1,:) = x0;
xdot_g = xdot0g;
% Initialize outputs
tp(1) = ti;  
xp(1,:) = x0; 
inputs_outputs = zeros(7,2);
N = length(x0);
Psi_int = eye(N);

%% Solve the ODEs
i = 1;
while tp(i) < tf
    hc = hlim;                % Reset current timestep to maximum allowable 
    relax_fac = 1;            % Reset relaxation factor
    eps = 10*RKFtol;          % Reset RKF epsilon
    b_it = 0;                 % Reset boundary iteration
    RKF_it = 0;               % Reset RKF approximations iteration
    while (eps > RKFtol || any(boundary < delta_b)) 
        if tc+hc > tf, hc = tf-tc; end  % Adjust for last time step
        xdot_it = 0;                    % Reset xdot iteration
        xdot_eps_rel = 10*xdotRelTol;   % Reset relative xdot epsilon
        xdot_eps_abs = 10*xdotAbsTol;   % Reset absolute xdot epsilon
        xdot_i = xdot_g;
        while xdot_eps_rel > xdotRelTol || xdot_eps_abs > xdotAbsTol
            [k1,inputs_outputs] = AEBL_dxdt_Sheng(tc,x(i,:),xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            x1 = x(i,:) + k1*hc/4;
            [k2,inputs_outputs] = AEBL_dxdt_Sheng(tc+hc/4,x1,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            x2 = x(i,:) + (3*k1/32+9/32*k2)*hc;
            [k3,inputs_outputs] = AEBL_dxdt_Sheng(tc+3/8*hc,x2,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            x3 = x(i,:) + (1932*k1-7200*k2+7296*k3)/2197*hc;
            [k4,inputs_outputs] = AEBL_dxdt_Sheng(tc+12/13*hc,x3,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            x4 = x(i,:) + (439/216*k1-8*k2+3680/513*k3-845/4104*k4)*hc;
            [k5,inputs_outputs] = AEBL_dxdt_Sheng(tc+hc,x4,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            x5 = x(i,:) + (-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)*hc;
            [k6,~] = AEBL_dxdt_Sheng(tc+hc/2,x5,xdot_g,tc,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            xdot_ng = (16/135*k1+6656/12825*k3+28561/56430*k4-9/50*k5+2/55*k6);
            x(i+1,:) = x(i,:) + xdot_ng*hc; 
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
        eps = sqrt((1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6)*(1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6)')*hc; % Error between 5th and 4th order approximations
        [boundary,inputs_outputs] = BL_boundaries(tc+hc,tc,x(i,9),x(i+1,9),x(i+1,10),x(i,10),x(i+1,12),x(i,12),inputs_outputs,mods,alpha1_0,c_n_alpha,TvL,q0,c_n1);
        % Check iterations on boundary
        if any(boundary < delta_b) 
            b_it = b_it+1;
            if b_it==b_max_it      
                break
            end
            hc = hc/2;
        % Check RKF iterations    
        elseif eps > RKFtol
            RKF_it = RKF_it+1;
            if RKF_it == RKF_it_max
                break
            end
            hc = hc/2;
        end
    end
    % Update calculation of Monodromy matrix
    Psi_int = monodromy_matrix_gen(Psi_method,Psi_int,N,tc,hc,x(i,:),xdot_i,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,xdotRelTol,xdot_max_it,dx);
    % Setup next time step
    tc = tc+hc;
    tp(i+1,1) = tc;      
    xp(i+1,:) = x(i+1,:);
    if i == 1
        xdot_0 = xdot_g;
    end
    i=i+1;
end

%% Truncate pre-allocated vectors
tp = tp(1:i); xp = xp(1:i,:);
