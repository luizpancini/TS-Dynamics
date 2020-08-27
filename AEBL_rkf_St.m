function [tp,xp,yp,xdotp] = AEBL_rkf_St(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,mods,q0,options)

%% Handle the inputs
if ~isempty(options)
    if any(strcmp(options,'hlim')) == 1, k=find(strcmp(options,'hlim')); hlim = options{k+1}; end
    if any(strcmp(options,'RKFtol')) == 1, k=find(strcmp(options,'RKFtol')); RKFtol = options{k+1}; end
    if any(strcmp(options,'xdotRelTol')) == 1, k=find(strcmp(options,'xdotRelTol')); xdotRelTol = options{k+1}; end
    if any(strcmp(options,'xdotAbsTol')) == 1, k=find(strcmp(options,'xdotAbsTol')); xdotAbsTol = options{k+1}; end
    if any(strcmp(options,'relax_fac_min')) == 1, k=find(strcmp(options,'relax_fac_min')); relax_fac_min = options{k+1}; end
    if any(strcmp(options,'LCO_tol')) == 1, k=find(strcmp(options,'LCO_tol')); LCO_tol = options{k+1}; end
    if any(strcmp(options,'SA_LCO_detect')) == 1, k=find(strcmp(options,'SA_LCO_detect')); SA_LCO_detect = options{k+1}; end
    if any(strcmp(options,'MA_LCO_detect')) == 1, k=find(strcmp(options,'MA_LCO_detect')); MA_LCO_detect = options{k+1}; end
    if any(strcmp(options,'min_SA_LCO_cycles')) == 1, k=find(strcmp(options,'min_SA_LCO_cycles')); min_SA_LCO_cycles = options{k+1}; end
    if any(strcmp(options,'min_MA_LCO_cycles')) == 1, k=find(strcmp(options,'min_MA_LCO_cycles')); min_MA_LCO_cycles = options{k+1}; end
    if any(strcmp(options,'alpha_bar_lim')) == 1, k=find(strcmp(options,'alpha_bar_lim')); alpha_bar_lim = options{k+1}; end
    if any(strcmp(options,'xdot_max_it')) == 1, k=find(strcmp(options,'xdot_max_it')); xdot_max_it = options{k+1}; end
    if any(strcmp(options,'b_max_it')) == 1, k=find(strcmp(options,'b_max_it')); b_max_it = options{k+1}; end
    if any(strcmp(options,'RKF_it_max')) == 1, k=find(strcmp(options,'RKF_it_max')); RKF_it_max = options{k+1}; end
end        

if ~exist('hlim','var'), hlim = 1e-4; end                               % Default maximum timestep
if ~exist('RKFtol','var'), RKFtol = 1e-8; end                          % Default RFK tolerance
if ~exist('xdotRelTol','var'), xdotRelTol = 1e-4; end                   % Default xdot relative tolerance
if ~exist('xdotAbsTol','var'), xdotAbsTol = 1e-4; end                   % Default xdot absolute tolerance
if ~exist('relax_fac_min','var'), relax_fac_min = 2^-3; end             % Default minimum relaxation factor (0-1): decrease for better (but slower) convergence on xdot
if ~exist('LCO_tol','var'), LCO_tol = 1e-7; end                         % Default steady-state LCO detection tolerance
if ~exist('SA_LCO_detect','var'), SA_LCO_detect = 1; end                % Default mode for single-amplitude LCO detection (0 = do not detect, 1 = detect)
if ~exist('MA_LCO_detect','var'), MA_LCO_detect = 1; end                % Default mode for multiple-amplitude LCO detection (0 = do not detect, 1 = detect)
if ~exist('min_SA_LCO_cycles','var'), min_SA_LCO_cycles = 20; end       % Default minimum single-amplitude LCO cycles to run
if ~exist('min_MA_LCO_cycles','var'), min_MA_LCO_cycles = 20; end       % Default minimum multiple-amplitude LCO cycles to run
if ~exist('alpha_bar_lim','var'), alpha_bar_lim = pi/4; end             % Default maximum effective AoA for DS model validity
if ~exist('xdot_max_it','var'), xdot_max_it = 100; end                  % Default maximum number of iterations for xdot implicit convergence
if ~exist('b_max_it','var'), b_max_it = 10; end                         % Default maximum number of iterations on boundary
if ~exist('RKF_it_max','var'), RKF_it_max = 10; end                     % Default maximum number of iterations for RKF approximations

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
xdot_g = xdot0;
% Initialize outputs
tp(1,1) = ti;  
xp(1,:) = x0; 
yp(1,:) = y0;
xdotp(1,:) = xdot0;    
% Initialize global variables
global continue_run divergent_flag inputs_outputs_save jM alpha_upeak alpha_lpeak h_upeak h_lpeak t_alpha_upeak t_alpha_lpeak t_h_upeak t_h_lpeak LCO_reached multipeaks_LCO it_multipeak_LCO 
if strcmp(continue_run,'n') || jM==1 
    inputs_outputs = zeros(7,2);
    it_multipeak_LCO = 0;
    divergent_flag = 0;
else
    inputs_outputs = inputs_outputs_save;
end
% Initialize LCO variables
LCO_alpha_u = 0; LCO_alpha_l = 0; LCO_xi_u = 0; LCO_xi_l = 0; it_upeak_alpha = 0; it_upeak_h = 0; it_lpeak_alpha = 0; it_lpeak_h = 0; 

%% Solve the ODEs
i = 1;
while tp(i) < tf
    hc = hlim;                % Reset current timestep to maximum allowable 
    relax_fac = 1;            % Reset relaxation factor
    b_it = 0;                 % Reset boundary iteration
    RKF_it = 0;               % Reset RKF approximations iteration
    eps = 10*RKFtol;          % Reset RKF epsilon
    while (eps > RKFtol || any(boundary < delta_b)) 
        if tc+hc > tf, hc = tf-tc; end  % Adjust for last time step
        xdot_it = 1;                    % Reset xdot iteration
        xdot_eps_rel = 10*xdotRelTol;   % Reset xdot epsilon
        xdot_eps_abs = 10*xdotAbsTol;   % Reset absolute xdot epsilon
        while xdot_eps_rel > xdotRelTol && xdot_eps_abs > xdotAbsTol 
            [k1,inputs_outputs] = AEBL_dxdt_standard(tc,x(i,:),xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof);
            x1 = x(i,:) + k1*hc/4;
            [k2,inputs_outputs] = AEBL_dxdt_standard(tc+hc/4,x1,xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof);
            x2 = x(i,:) + (3*k1/32+9/32*k2)*hc;
            [k3,inputs_outputs] =AEBL_dxdt_standard(tc+hc*3/8,x2,xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof);
            x3 = x(i,:) + (1932*k1-7200*k2+7296*k3)/2197*hc;
            [k4,inputs_outputs] = AEBL_dxdt_standard(tc+hc*12/13,x3,xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof);
            x4 = x(i,:) + (439/216*k1-8*k2+3680/513*k3-845/4104*k4)*hc;
            [k5,inputs_outputs] = AEBL_dxdt_standard(tc+hc,x4,xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof);
            x5 = x(i,:) + (-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)*hc;
            [k6,inputs_outputs] = AEBL_dxdt_standard(tc+hc/2,x5,xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof);
            xdot_ng = (16/135*k1+6656/12825*k3+28561/56430*k4-9/50*k5+2/55*k6);
            x(i+1,:) = x(i,:) + xdot_ng*hc; 
            [xdot_ip1g,inputs_outputs] = AEBL_dxdt_standard(tc+hc,x(i+1,:),xdot_ng,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof);
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
            if xdot_it==xdot_max_it, break, end
        end
        eps = sqrt((1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6)*(1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6)')*hc; % Error between 5th and 4th order approximations
        [boundary,inputs_outputs] = BL_boundaries(tc+hc,tc,x(i,9),x(i+1,10),x(i,10),x(i+1,12),x(i,12),inputs_outputs,mods,alpha1_0,c_n_alpha,TvL,q0,c_n1);
        if i>1
            boundary(11) = xdotp(i-1,15)*xdot_g(15);
            boundary(12) = xdotp(i-1,16)*xdot_g(16);
        end 
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
        if xdot_it==xdot_max_it
            warning(['xdotRelTol or xdotAbsTol not achieved at time ', num2str(tc,'%10.6f') ' s (xdot_eps_rel = ', num2str(xdot_eps_rel,'%.6f'), ', xdot_eps_abs = ', num2str(xdot_eps_abs,'%.6f'), ')'])
        end
 
        % Update derivatives
        xdotp(i,:) = xdot_g; % Update converged derivatives at current time step
        xdot_g = xdot_ip1g;  % Set initial guess for next time step
        % Check for divergent oscillations
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
        % Setup next time step
        inputs_outputs(3:end,1) = inputs_outputs(3:end,2);
        tc = tc+hc;
        tp(i+1,1) = tc;      
        xp(i+1,:) = x(i+1,:);
        yp(i+1,:) = AEBL_outputs_standard(tc,x(i+1,:),xdotp(i,:),inputs_outputs,U,b,ah,M,beta,c_n1,eta,E0,Df,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,b1,b2,b3,b4,A1,A2,A3,A4,c_m0,c_n_alpha,TvL,y0);
        i=i+1;
        if rem(i,100) == 0
            disp(['RKF45 progress: ',num2str((tc-ti)/(tf-ti)*100,'%10.2f') '%'])
        end
end
inputs_outputs_save = inputs_outputs;

% Truncate pre-allocated vectors
xdotp(i,:) = xdot_g;
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