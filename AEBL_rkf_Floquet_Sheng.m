function [tp,xp,xdotp,Psi_rkf,Psi_euler] = AEBL_rkf_Floquet_Sheng(tspan,x0,xdot0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options)
%% Handle the inputs
if ~isempty(options)
    if any(strcmp(options,'hlim')) == 1, k=find(strcmp(options,'hlim')); hlim = options{k+1}; end
    if any(strcmp(options,'RKFtol')) == 1, k=find(strcmp(options,'RKFtol')); RKFtol = options{k+1}; end
    if any(strcmp(options,'xdotRelTol')) == 1, k=find(strcmp(options,'xdotRelTol')); xdotRelTol = options{k+1}; end
    if any(strcmp(options,'xdotAbsTol')) == 1, k=find(strcmp(options,'xdotAbsTol')); xdotAbsTol = options{k+1}; end
    if any(strcmp(options,'relax_fac_min')) == 1, k=find(strcmp(options,'relax_fac_min')); relax_fac_min = options{k+1}; end
    if any(strcmp(options,'SS_abs_tol')) == 1, k=find(strcmp(options,'SS_abs_tol')); SS_abs_tol = options{k+1}; end
    if any(strcmp(options,'LCO_abs_tol')) == 1, k=find(strcmp(options,'LCO_abs_tol')); LCO_abs_tol = options{k+1}; end
    if any(strcmp(options,'LCO_rel_tol')) == 1, k=find(strcmp(options,'LCO_rel_tol')); LCO_rel_tol = options{k+1}; end
    if any(strcmp(options,'alpha_bar_lim')) == 1, k=find(strcmp(options,'alpha_bar_lim')); alpha_bar_lim = options{k+1}; end
    if any(strcmp(options,'xdot_max_it')) == 1, k=find(strcmp(options,'xdot_max_it')); xdot_max_it = options{k+1}; end
    if any(strcmp(options,'b_max_it')) == 1, k=find(strcmp(options,'b_max_it')); b_max_it = options{k+1}; end
    if any(strcmp(options,'RKF_it_max')) == 1, k=find(strcmp(options,'RKF_it_max')); RKF_it_max = options{k+1}; end
end        

if ~exist('hlim','var'), hlim = 1e-4; end                               % Default maximum timestep
if ~exist('RKFtol','var'), RKFtol = 1e-10; end                           % Default RFK tolerance
if ~exist('xdotRelTol','var'), xdotRelTol = 1e-12; end                   % Default xdot relative tolerance
if ~exist('xdotAbsTol','var'), xdotAbsTol = 1e-12; end                   % Default xdot absolute tolerance
if ~exist('relax_fac_min','var'), relax_fac_min = 2^-3; end             % Default minimum relaxation factor (0-1): decrease for better (but slower) convergence on xdot
if ~exist('SS_abs_tol','var'), SS_abs_tol = 1e-8; end                  % Default steady-state detection absolute tolerance
if ~exist('LCO_abs_tol','var'), LCO_abs_tol = 1e-6; end                 % Default single-amplitude LCO detection absolute tolerance
if ~exist('LCO_rel_tol','var'), LCO_rel_tol = 5e-4; end                 % Default single-amplitude LCO detection relative tolerance
if ~exist('alpha_bar_lim','var'), alpha_bar_lim = pi/4; end             % Default maximum effective AoA for DS model validity
if ~exist('xdot_max_it','var'), xdot_max_it = 5e1; end                  % Default maximum number of iterations for xdot implicit convergence
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
p_ratio = 1+sqrt(tf-ti);         % Pre-allocation ratio relative to constant time step - my simple formula
s0 = round(p_ratio*(tf-ti)/hlim);  % Initial size of the vectors
tp = zeros(s0,1); xp = zeros(s0,length(x0)); xdotp = zeros(s0,length(x0)); x = xp;
% Set initial conditions
x(1,:) = x0;
xdot_g = xdot0;
% Initialize outputs
tp(1,1) = ti;  
xp(1,:) = x0;
xdotp(1,:) = xdot0;
% Initialize global variables
global continue_run jT divergent_flag inputs_outputs_save jM LCO_reached LCO_period; LCO_reached = 0; 
% if ~strcmp(continue_run,'y') || jM==1 
%     inputs_outputs = zeros(7,2);
% else
    inputs_outputs = inputs_outputs_save;
% end
% LCO detection variables
N = size(x,2);
marker_state = 14;
min_crossings = 10;
x_max_reset = min_crossings;
x_at_marker_max = zeros(2*min_crossings,N); SS = 0; delta_x = 1; t_LCO_reached = 0; LCO_completed = 0; reset_ps = 1; x_max_count = 0; D_delta_x = 1;
% Floquet analysis 
Psi_rkf = eye(N);     % Initialize with identity matrix
Psi_euler = eye(N);   % Initialize with identity matrix
dx = 1e-12;           % Finite difference dx for Jacobian calculation

%% Solve the ODEs
i = 1;
n = 0;
while tp(i) < tf
    hc = hlim;                % Reset current timestep to maximum allowable   
    relax_fac = 1;            % Reset relaxation factor
    eps = 10*RKFtol;          % Reset RKF epsilon
    b_it = 0;                 % Reset boundary iteration
    RKF_it = 0;               % Reset RKF approximations iteration
    while (eps > RKFtol || any(boundary < delta_b)) 
        if tc+hc > tf, hc = tf-tc; end  % Adjust for last time step
        xdot_it = 1;                    % Reset xdot iteration
        xdot_eps_rel = 10*xdotRelTol;   % Reset xdot epsilon
        xdot_eps_abs = 10*xdotAbsTol;   % Reset absolute xdot epsilon
        inputs_outputs_s1 = inputs_outputs;
        while xdot_eps_rel > xdotRelTol || xdot_eps_abs > xdotAbsTol
            [k1,inputs_outputs] = AEBL_dxdt_Sheng(tc,x(i,:),xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            x1 = x(i,:) + k1*hc/4; inputs_outputs_s2 = inputs_outputs;
            [k2,inputs_outputs] = AEBL_dxdt_Sheng(tc+hc/4,x1,xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            x2 = x(i,:) + (3*k1/32+9/32*k2)*hc; inputs_outputs_s3 = inputs_outputs;
            [k3,inputs_outputs] = AEBL_dxdt_Sheng(tc+3/8*hc,x2,xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            x3 = x(i,:) + (1932*k1-7200*k2+7296*k3)/2197*hc; inputs_outputs_s4 = inputs_outputs;
            [k4,inputs_outputs] = AEBL_dxdt_Sheng(tc+12/13*hc,x3,xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            x4 = x(i,:) + (439/216*k1-8*k2+3680/513*k3-845/4104*k4)*hc; inputs_outputs_s5 = inputs_outputs;
            [k5,inputs_outputs] = AEBL_dxdt_Sheng(tc+hc,x4,xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            x5 = x(i,:) + (-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)*hc; inputs_outputs_s6 = inputs_outputs;
            [k6,inputs_outputs] = AEBL_dxdt_Sheng(tc+hc/2,x5,xdot_g,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            xdot_ng = (16/135*k1+6656/12825*k3+28561/56430*k4-9/50*k5+2/55*k6);
            x(i+1,:) = x(i,:) + xdot_ng*hc; 
            [xdot_ip1g,inputs_outputs] = AEBL_dxdt_Sheng(tc+hc,x(i+1,:),xdot_ng,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
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
            boundary(11) = xdotp(i-1,marker_state)*xdot_g(marker_state);
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
%             warning(['xdotRelTol or xdotAbsTol not achieved at time ', num2str(tc,'%10.6f') ' s (xdot_eps_rel = ', num2str(xdot_eps_rel), ', xdot_eps_abs = ', num2str(xdot_eps_abs), ')'])
        end
        % Update derivatives
        xdotp(i,:) = xdot_g; % Update converged derivatives at current time step
        xdot_g = xdot_ip1g;  % Set initial guess for next time step
        % Define a Poincaré section as when the marker state reaches a maximum 
        if i>1
        if xdotp(i-1,marker_state)>0 && xdotp(i,marker_state)<0 && jT>0 
            % Check for "steady-state" first
            if SS == 0
                n=n+1;
                x_at_marker_max(n,:) = x(i,:);
                if n == min_crossings
                    delta_x(1) = max(max(abs(diff(x_at_marker_max)))); % Maximum difference during the last n crossings
                elseif n == 2*min_crossings
                    delta_x(2) = max(max(abs(diff(x_at_marker_max(n-min_crossings+1:end,:))))); % Maximum difference during the last n crossings
                    D_delta_x_new = abs(diff(delta_x));
                    if D_delta_x_new < SS_abs_tol || abs(1-D_delta_x_new/D_delta_x) < 1e-3
                        SS = 1;
                        warning(['Steady-state started at time ', num2str(tc,'%.6f'),' s']);
                    else
                        n = min_crossings;
                        x_at_marker_max(1:n,:) = x_at_marker_max(n+1:2*min_crossings,:);
                        delta_x(1) = delta_x(2);
                        D_delta_x = D_delta_x_new;
                    end
                end
            % After "steady-state" is reached, set the Poincaré section
            else
                 x_max_count = x_max_count+1; if x_max_count==x_max_reset, x_max_count = 0; reset_ps = 1; end
                if reset_ps == 1
                    x_poincare = x(i,:);
                    reset_ps = 0;
                    LCO_reached = 0;
                    n = 0;
                    Psi_rkf = eye(N);     % Re-initialize with identity matrix
                    Psi_euler = eye(N);   % Re-initialize with identity matrix
                    warning(['Poincaré section defined at time ', num2str(tc,'%.6f'),' s']);
                elseif LCO_reached == 0 && max(abs(x(i,:)-x_poincare)) < LCO_abs_tol && max(abs(1-x(i,:)./x_poincare)) < LCO_rel_tol 
                    LCO_reached = 1;
                    t_LCO_reached = tc;
                    warning(['LCO started at time ', num2str(tc,'%.6f'),' s']);
                end
                % Marker state reaches a maximum while piercing the Poincaré
                % section at the same point
                if LCO_reached == 1 && tc > t_LCO_reached && max(abs(x(i,:)-x_poincare)) < LCO_abs_tol && max(abs(1-x(i,:)./x_poincare)) < LCO_rel_tol %&& max(abs(xdotp(i+1,:)-xdot_poincare)) < LCO_abs_tol 
                    LCO_completed = 1;
                    warning(['LCO finished at time ', num2str(tc,'%.6f'),' s']);
                    LCO_period = tc - t_LCO_reached;             
                end
            end
        end
        end
        % Floquet analysis - if the LCO is in its second cycle after it has been detected, calculate Jacobian and increment Monodromy matrix
        if LCO_reached == 1 
            % Apply Euler method to calculate the Monodromy matrix
            J = zeros(N);
            for j=1:N
                Dx = zeros(1,N); Dx(j) = dx;
                [F_pDx,~] = AEBL_dxdt_Sheng(tc,x(i,:)+Dx,xdotp(i,:),inputs_outputs_s1,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                [F_mDx,~] = AEBL_dxdt_Sheng(tc,x(i,:)-Dx,xdotp(i,:),inputs_outputs_s1,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                J(:,j) = (F_pDx-F_mDx)/(2*dx);
            end
            Psi_euler = Psi_euler+hc*J*Psi_euler;
            % Apply RKF to calculate the Monodromy matrix
            % Step 1 - Jacobian is the same from Euler method for step 1
            k1_J = J*Psi_rkf;
            Psi_bar_1 = Psi_rkf+k1_J*hc/4;
            % Step 2
            J = zeros(N);
            for j=1:N
                Dx = zeros(1,N); Dx(j) = dx;
                [F_pDx,~] = AEBL_dxdt_Sheng(tc+hc/4,x1+Dx,xdotp(i,:),inputs_outputs_s2,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                [F_mDx,~] = AEBL_dxdt_Sheng(tc+hc/4,x1-Dx,xdotp(i,:),inputs_outputs_s2,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                J(:,j) = (F_pDx-F_mDx)/(2*dx);
            end
            k2_J = J*Psi_bar_1;
            Psi_bar_2 = Psi_rkf + (3*k1_J/32+9/32*k2_J)*hc;
            % Step 3
            J = zeros(N);
            for j=1:N
                Dx = zeros(1,N); Dx(j) = dx;
                [F_pDx,~] = AEBL_dxdt_Sheng(tc+3/8*hc,x2+Dx,xdotp(i,:),inputs_outputs_s3,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                [F_mDx,~] = AEBL_dxdt_Sheng(tc+3/8*hc,x2-Dx,xdotp(i,:),inputs_outputs_s3,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                J(:,j) = (F_pDx-F_mDx)/(2*dx);
            end
            k3_J = J*Psi_bar_2;
            Psi_bar_3 = Psi_rkf + (1932*k1_J-7200*k2_J+7296*k3_J)/2197*hc;
            % Step 4
            J = zeros(N);
            for j=1:N
                Dx = zeros(1,N); Dx(j) = dx;
                [F_pDx,~] = AEBL_dxdt_Sheng(tc+12/13*hc,x3+Dx,xdotp(i,:),inputs_outputs_s4,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                [F_mDx,~] = AEBL_dxdt_Sheng(tc+12/13*hc,x3-Dx,xdotp(i,:),inputs_outputs_s4,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                J(:,j) = (F_pDx-F_mDx)/(2*dx);
            end
            k4_J = J*Psi_bar_3;
            Psi_bar_4 = Psi_rkf + (439/216*k1_J-8*k2_J+3680/513*k3_J-845/4104*k4_J)*hc;
            % Step 5
            J = zeros(N);
            for j=1:N
                Dx = zeros(1,N); Dx(j) = dx;
                [F_pDx,~] = AEBL_dxdt_Sheng(tc+hc,x4+Dx,xdotp(i,:),inputs_outputs_s5,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                [F_mDx,~] = AEBL_dxdt_Sheng(tc+hc,x4-Dx,xdotp(i,:),inputs_outputs_s5,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                J(:,j) = (F_pDx-F_mDx)/(2*dx);
            end
            k5_J = J*Psi_bar_4;
            Psi_bar_5 = Psi_rkf + (-8/27*k1_J+2*k2_J-3544/2565*k3_J+1859/4104*k4_J-11/40*k5_J)*hc;
            % Step 6
            J = zeros(N);
            for j=1:N
                Dx = zeros(1,N); Dx(j) = dx;
                [F_pDx,~] = AEBL_dxdt_Sheng(tc+hc/2,x5+Dx,xdotp(i,:),inputs_outputs_s6,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                [F_mDx,~] = AEBL_dxdt_Sheng(tc+hc/2,x5-Dx,xdotp(i,:),inputs_outputs_s6,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
                J(:,j) = (F_pDx-F_mDx)/(2*dx);
            end
            k6_J = J*Psi_bar_5;
            Psi_rkf = Psi_rkf + (16/135*k1_J+6656/12825*k3_J+28561/56430*k4_J-9/50*k5_J+2/55*k6_J)*hc;
        end
        % Setup next time step
        inputs_outputs(3:end,1) = inputs_outputs(3:end,2);
        tc = tc+hc;
        tp(i+1,1) = tc;      
        xp(i+1,:) = x(i+1,:);
        i=i+1;
        if rem(i,1e4) == 0
            disp(['RKF45 progress: ',num2str((tc-ti)/(tf-ti)*100,'%10.2f') '%'])
        end
        if abs(inputs_outputs(4,2)) > alpha_bar_lim
            warning(['Effective AoA > ', num2str(alpha_bar_lim*180/pi,'%10.2f'),'° - stopping integration'])
            divergent_flag = 1;
            break
        end   
        if LCO_completed == 1
            break
        end
end
inputs_outputs_save = inputs_outputs;

% Truncate pre-allocated vectors
xdotp(i,:) = xdot_g;
tp = tp(1:i); xp = xp(1:i,:); xdotp = xdotp(1:i,:);

% figure;plot(tp,xp(:,14));grid
end