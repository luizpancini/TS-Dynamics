function [tp,xp,yp,xdotp] = AEBL_rkf_Sheng_no_accel_terms(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options)

%% Handle the inputs
if ~isempty(options)
    if any(strcmp(options,'hlim')) == 1, k=find(strcmp(options,'hlim')); hlim = options{k+1}; end
    if any(strcmp(options,'RKFtol')) == 1, k=find(strcmp(options,'RKFtol')); RKFtol = options{k+1}; end
    if any(strcmp(options,'LCO_tol')) == 1, k=find(strcmp(options,'LCO_tol')); LCO_tol = options{k+1}; end
    if any(strcmp(options,'SA_LCO_detect')) == 1, k=find(strcmp(options,'SA_LCO_detect')); SA_LCO_detect = options{k+1}; end
    if any(strcmp(options,'MA_LCO_detect')) == 1, k=find(strcmp(options,'MA_LCO_detect')); MA_LCO_detect = options{k+1}; end
    if any(strcmp(options,'min_SA_LCO_cycles')) == 1, k=find(strcmp(options,'min_SA_LCO_cycles')); min_SA_LCO_cycles = options{k+1}; end
    if any(strcmp(options,'min_MA_LCO_cycles')) == 1, k=find(strcmp(options,'min_MA_LCO_cycles')); min_MA_LCO_cycles = options{k+1}; end
    if any(strcmp(options,'alpha_bar_lim')) == 1, k=find(strcmp(options,'alpha_bar_lim')); alpha_bar_lim = options{k+1}; end
    if any(strcmp(options,'b_max_it')) == 1, k=find(strcmp(options,'b_max_it')); b_max_it = options{k+1}; end
    if any(strcmp(options,'RKF_it_max')) == 1, k=find(strcmp(options,'RKF_it_max')); RKF_it_max = options{k+1}; end
end        

if ~exist('hlim','var'), hlim = 1e-4; end                               % Default maximum timestep
if ~exist('RKFtol','var'), RKFtol = 1e-8; end                          % Default RFK tolerance
if ~exist('LCO_tol','var'), LCO_tol = 1e-7; end                         % Default steady-state LCO detection tolerance
if ~exist('SA_LCO_detect','var'), SA_LCO_detect = 1; end                % Default mode for single-amplitude LCO detection (0 = do not detect, 1 = detect)
if ~exist('MA_LCO_detect','var'), MA_LCO_detect = 1; end                % Default mode for multiple-amplitude LCO detection (0 = do not detect, 1 = detect)
if ~exist('min_SA_LCO_cycles','var'), min_SA_LCO_cycles = 20; end       % Default minimum single-amplitude LCO cycles to run
if ~exist('min_MA_LCO_cycles','var'), min_MA_LCO_cycles = 20; end       % Default minimum multiple-amplitude LCO cycles to run
if ~exist('alpha_bar_lim','var'), alpha_bar_lim = pi/4; end             % Default maximum effective AoA for DS model validity
if ~exist('b_max_it','var'), b_max_it = 12; end                         % Default maximum number of iterations on boundary
if ~exist('RKF_it_max','var'), RKF_it_max = 12; end                     % Default maximum number of iterations for RKF approximations

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
tp = zeros(s0,1); xp = zeros(s0,length(x0)); yp = zeros(s0,length(y0)); xdotp = zeros(s0,length(x0)); x = xp;
% Set initial conditions
x(1,:) = x0;
xdot = xdot0;
% Initialize outputs
tp(1,1) = ti;  
xp(1,:) = x0; 
yp(1,:) = y0;
xdotp(1,:) = xdot0;    
% Initialize global variables
global continue_run divergent_flag inputs_outputs_save jM alpha_upeak alpha_lpeak h_upeak h_lpeak t_alpha_upeak t_alpha_lpeak t_h_upeak t_h_lpeak LCO_reached multipeaks_LCO it_multipeak_LCO LCO_cycles
if ~strcmp(continue_run,'y') || jM==1 
    inputs_outputs = zeros(7,2);
else
    inputs_outputs = inputs_outputs_save;
end
% Initialize LCO variables
LCO_alpha_u = 0; LCO_alpha_l = 0; LCO_xi_u = 0; LCO_xi_l = 0; it_upeak_alpha = 0; it_upeak_h = 0; it_lpeak_alpha = 0; it_lpeak_h = 0; LCO_cycles = min_SA_LCO_cycles;

%% Solve the ODEs
i = 1;
while tp(i) < tf
    hc = hlim;                % Reset current timestep to maximum allowable 
    b_it = 0;                 % Reset boundary iteration
    RKF_it = 0;               % Reset RKF approximations iteration
    while (eps > RKFtol || any(boundary < delta_b)) 
        if tc+hc > tf, hc = tf-tc; end  % Adjust for last time step
        [k1,inputs_outputs] = AEBL_dxdt_Sheng_no_accel_terms(tc,x(i,:),inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x1 = x(i,:) + k1*hc/4;
        [k2,inputs_outputs] = AEBL_dxdt_Sheng_no_accel_terms(tc+hc/4,x1,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x2 = x(i,:) + (3*k1/32+9/32*k2)*hc;
        [k3,inputs_outputs] = AEBL_dxdt_Sheng_no_accel_terms(tc+3/8*hc,x2,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x3 = x(i,:) + (1932*k1-7200*k2+7296*k3)/2197*hc;
        [k4,inputs_outputs] = AEBL_dxdt_Sheng_no_accel_terms(tc+12/13*hc,x3,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x4 = x(i,:) + (439/216*k1-8*k2+3680/513*k3-845/4104*k4)*hc;
        [k5,inputs_outputs] = AEBL_dxdt_Sheng_no_accel_terms(tc+hc,x4,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x5 = x(i,:) + (-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)*hc;
        [k6,inputs_outputs] = AEBL_dxdt_Sheng_no_accel_terms(tc+hc/2,x5,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        x(i+1,:) = x(i,:) + (16/135*k1+6656/12825*k3+28561/56430*k4-9/50*k5+2/55*k6)*hc;
        [xdot,inputs_outputs] = AEBL_dxdt_Sheng_no_accel_terms(tc+hc,x(i+1,:),inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
        eps = sqrt((1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6)*(1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6)')*hc; % Error between 5th and 4th order approximations
        [boundary,inputs_outputs] = BL_boundaries(tc+hc,tc,x(i,9),x(i+1,10),x(i,10),x(i+1,12),x(i,12),inputs_outputs,mods,alpha1_0,c_n_alpha,TvL,q0,c_n1);
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
        % Setup next time step
        inputs_outputs(3:end,1) = inputs_outputs(3:end,2);
        tc = tc+hc;
        tp(i+1,1) = tc;      
        xp(i+1,:) = x(i+1,:);
        yp(i+1,:) = AEBL_outputs_Sheng_no_accel_terms(tc,x(i+1,:),inputs_outputs,U,b,ah,M,beta,eta,E0,Df,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,b1,b2,b3,b4,A1,A2,A3,A4,c_m0,c_n_alpha,TvL,alpha_ss,alpha_ds0,q0,alpha1_0,S1,S2,y0);
        xdotp(i+1,:) = xdot;
        i=i+1;
        if rem(i,100) == 0
            disp(['RKF45 progress: ',num2str(tc/tf*100,'%10.2f') '%'])
        end
        if abs(inputs_outputs(4,2)) > alpha_bar_lim
            warning(['Effective AoA > ', num2str(alpha_bar_lim*180/pi,'%10.2f'),'° - stopping integration'])
            divergent_flag = 1;
            break
        end
        % Peaks finder
        if SA_LCO_detect == 1 && i>1 % Check if has reached LCOs
            if xdotp(i-1,16) > 0 && xdotp(i,16) < 0 % alpha reaches a maximum
                it_upeak_alpha = it_upeak_alpha+1;
                alpha_upeak(it_upeak_alpha) = x(i,16);
                t_alpha_upeak(it_upeak_alpha) = tc;
                if it_upeak_alpha > min_SA_LCO_cycles
                    if max(abs(diff(alpha_upeak(end-min_SA_LCO_cycles+1:end))))<LCO_tol
                        LCO_alpha_u = 1;
                    end
                end
                if MA_LCO_detect == 1
                    if it_upeak_alpha > min_MA_LCO_cycles && alpha_upeak(it_upeak_alpha)<alpha_upeak(it_upeak_alpha-1) && alpha_upeak(it_upeak_alpha-1)>alpha_upeak(it_upeak_alpha-2)
                          it_multipeak_LCO = it_multipeak_LCO+1;
                    end
                end
            elseif xdotp(i-1,16) < 0 && xdotp(i,16) > 0 % alpha reaches a minimum
                it_lpeak_alpha = it_lpeak_alpha+1;
                alpha_lpeak(it_lpeak_alpha) = x(i,16);
                t_alpha_lpeak(it_lpeak_alpha) = tc;
                if it_lpeak_alpha > min_SA_LCO_cycles
                    if max(abs(diff(alpha_lpeak(end-min_SA_LCO_cycles+1:end))))<LCO_tol
                        LCO_alpha_l = 1;
                    end
                end
            end
            if xdotp(i-1,15) > 0 && xdotp(i,15) < 0 % h reaches a maximum
                it_upeak_h = it_upeak_h+1;
                h_upeak(it_upeak_h) = x(i,15);
                t_h_upeak(it_upeak_h) = tc;
                if it_upeak_h > min_SA_LCO_cycles
                    if max(abs(diff(h_upeak(end-min_SA_LCO_cycles+1:end))))<LCO_tol
                        LCO_xi_u = 1; 
                    end
                end
            elseif xdotp(i-1,15) < 0 && xdotp(i,15) > 0 % h reaches a minimum
                it_lpeak_h = it_lpeak_h+1;
                h_lpeak(it_lpeak_h) = x(i,15);
                t_h_lpeak(it_lpeak_h) = tc;
                if it_lpeak_h > min_SA_LCO_cycles
                    if max(abs(diff(h_lpeak(end-min_SA_LCO_cycles+1:end))))<LCO_tol
                        LCO_xi_l = 1; 
                    end
                end
            end
            if LCO_alpha_u == 1 && LCO_alpha_l == 1 && LCO_xi_u == 1 && LCO_xi_l == 1 && isempty(LCO_reached)
                LCO_reached = 1;
                warning('Steady-state LCOs detected - stopping at the end of this run');
            end
            if MA_LCO_detect == 1 && it_multipeak_LCO == min_MA_LCO_cycles && isempty(LCO_reached)       
                multipeaks_LCO = 1;
                LCO_reached = 1;
                warning(['Multipeak LCOs detected and at least ', num2str(min_MA_LCO_cycles),' cycles completed - stopping at the end of this run'])
            end
        end
end
inputs_outputs_save = inputs_outputs;

% Truncate pre-allocated vectors
tp = tp(1:i); xp = xp(1:i,:); yp = yp(1:i,:); xdotp = xdotp(1:i,:);

% Truncate global vectors for only last min_SA_LCO_cycles if the LCOs are
% single-amplitude
if ~isempty(LCO_reached) && isempty(multipeaks_LCO)
    alpha_upeak = alpha_upeak(end-min_SA_LCO_cycles+1:end);
    t_alpha_upeak = t_alpha_upeak(end-min_SA_LCO_cycles+1:end);
    alpha_lpeak = alpha_lpeak(end-min_SA_LCO_cycles+1:end);
    t_alpha_lpeak = t_alpha_lpeak(end-min_SA_LCO_cycles+1:end);
    h_upeak = h_upeak(end-min_SA_LCO_cycles+1:end);
    t_h_upeak = t_h_upeak(end-min_SA_LCO_cycles+1:end);
    h_lpeak = h_lpeak(end-min_SA_LCO_cycles+1:end);
    t_h_lpeak = t_h_lpeak(end-min_SA_LCO_cycles+1:end);
end