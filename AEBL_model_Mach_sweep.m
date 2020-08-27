clc
clear all
close all

%% Inputs
h = 0; [a_inf,rho,~,~] = std_atmosphere_calc(h); a_inf = 343; rho = 1.225;
typsec = 'mine';
airfoil = 'NACA0012';

%% TS flutter data
[M_F,omega_F] = flutter_data(a_inf,rho,typsec,airfoil);
if strcmp(typsec,'Fung')
    M_F = 0.07693; % Override (already been using this value)
end

%% Mach sweep options
Mach_init = 0.10952;        % Starting M
dM = -1e-5;                  % M step
rundown = 0;                % Choose to restart going in the opposite direction after divergent osc. are reached
Mach_limit = 'n';           % Set a limit Mach number to stop
Mach_final = 0.07750;       % Value of Mach_limit
multipeak_LCO_tspan = 30;   % Time to run after multiple-amplitude LCOs are detected
singlepeak_LCO_tspan = 5;  % Time to run after single-amplitude LCOs are detected
tf0 = 10;                   % Minimum time of first run (to discard as transients)
ts = 5;                     % Maximum single runtime 
load_inicon_data = 'y';     % Option to load initial conditions from file

%% Run the Mach sweep
LCO_flag = 0; fixed_point_flag = 0; Mach_limit_reached = 0;
global jM divergent_flag it_multipeak_LCO continue_run; jM = 0; it_multipeak_LCO = 0; divergent_flag = 0; continue_run = 'y';
global inputs_outputs_save alpha_upeak alpha_lpeak h_upeak h_lpeak t_alpha_upeak t_alpha_lpeak t_h_upeak t_h_lpeak LCO_reached multipeaks_LCO

while fixed_point_flag == 0
    if divergent_flag == 1 && rundown == 1 % Upper limit has been reached - begin decreasing airspeed sweep from Mach_init
        LCO_flag = 0; jM = 0; it_multipeak_LCO = 0; divergent_flag = 0;
        dM = -dM;
    elseif (divergent_flag == 1 && rundown == 0) || Mach_limit_reached == 1
        fixed_point_flag = 1; % Just to stop the program
    end
    M = Mach_init;
while divergent_flag == 0 && fixed_point_flag == 0
    if LCO_flag == 1 % If LCO has been reached, advance the Mach
        M = M+dM;
    end
    M
    tf = tf0;
    tf_save = 0;
    final_run = 0;

    while tf_save < tf
        %% Reset total vectors
        t_total = [];
        alpha_total = []; 
        h_b_total = [];
        alpha_upeak_total = [];
        alpha_lpeak_total = [];
        t_alpha_upeak_total = [];
        t_alpha_lpeak_total = [];
        h_upeak_total = [];
        h_lpeak_total = [];
        t_h_upeak_total = [];
        t_h_lpeak_total = [];

        %% Loop until LCOs are detected or for at least the discard time
        while isempty(LCO_reached) || tf_save < tf
            
            jM = jM+1;
            
            %% Test conditions and flow
            dt = 5e-5;
            U = M*a_inf;    
            qinf = 1/2*rho*U^2;
            beta = sqrt(1-M^2);

            %% Choose if modifications will be allowed - 'on' for Sheng's modification, 'off' for standard BL
            mods = 'on'; 

            %% Airfoil structural parameters
            [b,ah,m,I_alpha,S_alpha,Kh,Ka,Ch,Ca,Kh3,Ka3,h_I,alpha_I,plunge_dof] = typsec_parameters(typsec,rho);
            typsec_filename = get_typsec_filename(typsec,b,ah,m,I_alpha,S_alpha,Kh,Ka,Ka3,h_I,alpha_I,rho,h);

            %% Airfoil aerodynamic parameters 
            [c_n_alpha,S1,S2,alpha1_0,alpha_ss,delta_alpha1,alpha_ds0,Df,TvL,c_n1,K0,K1,K2,x_ac,kappa,Tp,Ta,Tf0,Tv0,q0,E0,c_d0,c_m0,eta] = airfoil_parameters(airfoil,M,U,b);

            %% Indicial response parameters
            % Circulatory and non-circulatory unit step AoA indicial response - BL model
            [A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I] = read_indicial_params(M,beta,b,a_inf);

            %% State space matrices
            [A,B] = get_SS_matrices_AEBL(U,b,beta,M,a_inf);

            %% Initial conditions
            if strcmp(continue_run,'y') && jM>1
                x0 = x0_save;
                xdot0 = xdot0_save;
                tspan = [tf_save tf_save+ts]
            elseif strcmp(load_inicon_data,'y')
                filename = [typsec_filename,'_M_', num2str(M,'%.5f'), '.mat']; 
                load(filename,'x0_save','xdot0_save','y0_save','inputs_outputs_save')
                x0 = x0_save;
                xdot0 = xdot0_save;
                y0 = y0_save;
                tspan = [0 ts]
            else
                % Specified initial conditions
                h_0 = 0.00*b; 
                h_dot0 = 0; 
                h_2dot0 = 0; 
                alpha_0 = 1*pi/180; 
                alpha_dot0 = 0; 
                alpha_2dot0 = 0;
                % Initialize BL and structural states 
                [x0,y0,xdot0] = initialize_BL(h_0,h_dot0,h_2dot0,alpha_0,alpha_dot0,alpha_2dot0,A,B,U,M,b,ah,alpha1_0,delta_alpha1,S1,S2,c_n1,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,eta,E0,Df,alpha_ds0,alpha_ss,q0,mods);
                tspan = [0 ts]
            end

            %% ODE solver
            RKF_relax_fac_min = 1/2^0; RKFtol = 1e-10; xdotRelTol = 1e-12; xdotAbsTol = 1e-12; xdot_max_it = 10; disp_boundaries = 1; alpha_bar_lim = 60*pi/180;
            options = {'SA_LCO_detect',1,'MA_LCO_detect',1,'hlim',dt,'relax_fac_min',RKF_relax_fac_min,'RKFtol',RKFtol,'xdotRelTol',xdotRelTol,'xdotAbsTol',xdotAbsTol,'xdot_max_it',xdot_max_it,'disp_boundaries',disp_boundaries,'alpha_bar_lim',alpha_bar_lim};
            if strcmp(mods,'off') % Standard BL
                [t,x,y,xdot] = AEBL_rkf_St(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_I,Ka,Ka3,alpha_I,Ch,Ca,plunge_dof,mods,q0,options);
            else                  % Sheng's modifications
                [t,x,y,xdot] = AEBL_rkf_Sheng(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_I,Ka,Ka3,alpha_I,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
            end

            %% Save final conditions for next run's initial conditions
            x0_save = x(end,:);
            y0_save = y(end,:);
            xdot0_save = xdot(end,:);
            tf_save = t(end);

            %% Output variables 
            h_b = x(:,15)/b;
            alpha = x(:,16);

            %% Increment total vectors
            t_total = [t_total; t];
            alpha_total = [alpha_total; alpha]; 
            h_b_total = [h_b_total; h_b];
            alpha_upeak_total = [alpha_upeak_total; alpha_upeak'];
            alpha_lpeak_total = [alpha_lpeak_total; alpha_lpeak'];
            t_alpha_upeak_total = [t_alpha_upeak_total; t_alpha_upeak'];
            t_alpha_lpeak_total = [t_alpha_lpeak_total; t_alpha_lpeak'];
            h_upeak_total = [h_upeak_total; h_upeak'];
            h_lpeak_total = [h_lpeak_total; h_lpeak'];
            t_h_upeak_total = [t_h_upeak_total; t_h_upeak'];
            t_h_lpeak_total = [t_h_lpeak_total; t_h_lpeak'];

            %% Check for fixed point
            alpha_max = max(abs(alpha));
            if abs(alpha_I-alpha_max) < 1e-2
                warning('Oscillations too close to wind-off angle reached - assuming fixed point has been reached and stopping')
                fixed_point_flag = 1;
                break;
            end

            %% Check divergent oscillations
            if divergent_flag == 1
                break
            end
            
            %% Prevent short time span from not detecting MA LCOs 
            if it_multipeak_LCO >= 50
                LCO_reached = 1;
                multipeaks_LCO = 1;
            end
        end

        %% Check divergent oscillations
        if divergent_flag == 1
            break
        end

        %% Check if this was the final run
        if final_run == 1
            break
        end

        %% Check for LCOs type and set final time span accordingly
        if it_multipeak_LCO>=10 
            warning(['Multipeak LCOs have been detected - continuing from last run and setting time span to ', num2str(multipeak_LCO_tspan,'%.f') ' seconds'] )
            tf = multipeak_LCO_tspan+tf_save; 
            final_run = 1;
            alpha_upeak = []; alpha_lpeak = []; h_upeak = []; h_lpeak = []; t_alpha_upeak = []; t_alpha_lpeak = []; t_h_upeak = []; t_h_lpeak = []; LCO_reached  = []; multipeaks_LCO = 0; it_multipeak_LCO = 0;
        elseif it_multipeak_LCO<10
            warning(['Single-amplitude LCOs have been detected - continuing from last run and setting time span to ', num2str(singlepeak_LCO_tspan,'%.f') ' seconds'] )
            tf = singlepeak_LCO_tspan+tf_save; 
            final_run = 1;
            alpha_upeak = []; alpha_lpeak = []; h_upeak = []; h_lpeak = []; t_alpha_upeak = []; t_alpha_lpeak = []; t_h_upeak = []; t_h_lpeak = []; LCO_reached  = []; multipeaks_LCO = 0; it_multipeak_LCO = 0;
        end

    end

    %% Save LCO data
    % Rename to total vectors
    t = t_total;
    alpha = alpha_total; 
    h_b = h_b_total;
    alpha_upeak = alpha_upeak_total;
    alpha_lpeak = alpha_lpeak_total;
    t_alpha_upeak = t_alpha_upeak_total;
    t_alpha_lpeak = t_alpha_lpeak_total;
    h_upeak = h_upeak_total;
    h_lpeak = h_lpeak_total;
    t_h_upeak = t_h_upeak_total;
    t_h_lpeak = t_h_lpeak_total;
    % Save    
    if divergent_flag == 0 
        LCO_flag = 1; 
        if multipeaks_LCO == 0 % For single-amplitude LCOs
            % Get the reduced frequency and phase angle of motion
            [k,phi] = get_LCO_k_and_phi(b,U,alpha,h_b,h_upeak/b,h_lpeak/b,t_alpha_upeak,t_alpha_lpeak,t_h_upeak,t_h_lpeak);
            % Save only last value as LCO amplitude
            alpha_LCO_u = alpha_upeak(end); t_alpha_LCO_u = t_alpha_upeak(end);
            alpha_LCO_l = alpha_lpeak(end); t_alpha_LCO_l = t_alpha_lpeak(end);
            hb_LCO_u = h_upeak(end)/b; t_hb_LCO_u = t_h_upeak(end);
            hb_LCO_l = h_lpeak(end)/b; t_hb_LCO_l = t_h_lpeak(end);
        else                       % For multiple-amplitude LCOs
            alpha_LCO_u = alpha_upeak; t_alpha_LCO_u = t_alpha_upeak;
            alpha_LCO_l = alpha_lpeak; t_alpha_LCO_l = t_alpha_lpeak;
            hb_LCO_u = h_upeak/b; t_hb_LCO_u = t_h_upeak;
            hb_LCO_l = h_lpeak/b; t_hb_LCO_l = t_h_lpeak;
            k = nan; phi = nan;
        end
        % Save output variables 
        filename = [typsec_filename,'_M_', num2str(M,'%.5f'), '_upper_branch.mat']; 
        typsec_filepath = filename;
        save (typsec_filepath,'x0_save','xdot0_save','tf_save','y0_save','inputs_outputs_save','alpha_LCO_u','alpha_LCO_l','hb_LCO_u','hb_LCO_l','t_alpha_LCO_u','t_alpha_LCO_l','t_hb_LCO_u','t_hb_LCO_l','M','multipeaks_LCO','M_F','k','phi')
        it_multipeak_LCO = 0;
        % Plots
        set(0,'DefaultTextInterpreter','latex')
        set(0,'DefaultLegendInterpreter','latex')
        axes_size = 20;
        lw = 0.75;
        close all   
        figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
        axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman');
        hold(axes1,'on'); 
        plot(t,alpha*180/pi,'k-',t_alpha_LCO_u,alpha_LCO_u*180/pi,'ko',t_alpha_LCO_l,alpha_LCO_l*180/pi,'ko','LineWidth',lw,'Parent',axes1);
        ylabel('$\alpha$ [deg]','FontWeight','normal','FontSize',axes_size);
        xlabel('t [s]','FontWeight','normal','FontSize',axes_size);
        grid on
        title(['M = ', num2str(M)],'FontWeight','normal','FontSize',axes_size);
        drawnow 
    else
        LCO_flag = 0;
    end
    alpha_upeak = []; alpha_lpeak = []; h_upeak = []; h_lpeak = []; t_alpha_upeak = []; t_alpha_lpeak = []; t_h_upeak = []; t_h_lpeak = []; LCO_reached  = []; multipeaks_LCO = 0;

    %% Limit Mach number
    if strcmp(Mach_limit,'y')
        if M+dM >= Mach_final
            Mach_limit_reached = 1;
            break
        end
    end

end
end