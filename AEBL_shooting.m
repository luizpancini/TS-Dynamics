clc
clear all
close all

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
axes_size = 16;
lw = 0.75;

M_F = 0.07693;
M_init = 0.07692;
M_final = 0.07691;
dM = -1e-5;
M_vec = M_init:dM:M_final;
load_from_mat = 0; % Choose to load from (1, 2 or 3) available saved data or (0) manually specified initial conditions
global inputs_outputs_save

for i=1:length(M_vec)
    M = M_vec(i)
    
    %% Initial conditions
    if i == 1 
        filename1 = ['typsec_Fung_M_',num2str(M,'%.5f'),'_LCO_data.mat'];
        filename2 = ['typsec_Fung_M_',num2str(M,'%.6f'),'_Shooting_data_new.mat'];
        if load_from_mat == 1
            load(filename1, 'x0_save','xdot0_save','inputs_outputs_save','LCO_period');
            x0 = x0_save;
            xdot0 = xdot0_save;
            T0 = LCO_period;
        elseif load_from_mat == 2
            load(filename2,'x0','xdot0','inputs_outputs_save','T0');
        else
            x0 = [-0.000290753299112  -0.000153514740679  -0.000010190063473   0.000002515862408  -0.000002539551493  -0.000001013312792 -0.000002286249264   0.000001235933812  -0.012280449377022   0.998078591998070                   0   0.996031560793824 0.016916023556038   0.537802222909707  -0.000182462636501  -0.016548504235627];
            xdot0 = [-0.006130871522298   0.002265493664308   0.000263321915156   0.000301424131384   0.000069941668242   0.000028261173736 0.004788275811627   0.000141599779102  -0.241454077226535  -0.000666745579136                   0   0.007587545519224 0.565281721397623  51.274342717403776   0.016916029072635   0.537802723635741];
            T0 = 0.112877778431930;   
            inputs_outputs_save = zeros(7,2);
        end
    else % Take last converged values as initial conditions 
        filename2 = ['typsec_Fung_M_',num2str(M,'%.6f'),'_Shooting_data_new.mat'];
        if load_from_mat == 2 && exist(filename2,'file')
            load(filename2,'x0','xdot0','inputs_outputs_save','T0');
        end
    end


    %% Test conditions and airfoil data
    dt = 1e-4;
    h = 0; [a_inf,rho,~,~] = std_atmosphere_calc(h); a_inf = 343; rho = 1.225;
    U = M*a_inf;    
    qinf = 1/2*rho*U^2;
    beta = sqrt(1-M^2);

    % Choose if modifications will be allowed - 'on' for Sheng's modification, 'off' for standard BL
    mods = 'on'; 

    % Choose integration options
    RKF_relax_fac_min = 1/2^0; RKFtol = 1e-10; xdotRelTol = 1e-6; xdotAbsTol = 1e-6; xdot_max_it = 10;
    options = {'hlim',dt,'relax_fac_min',RKF_relax_fac_min,'RKFtol',RKFtol,'xdotRelTol',xdotRelTol,'xdotAbsTol',xdotAbsTol,'xdot_max_it',xdot_max_it};

    % Airfoil structural parameters
    typsec = 'Fung';
    [b,ah,m,I_alpha,S_alpha,Kh,Ka,Ch,Ca,Kh3,Ka3,h_F,alpha_F,plunge_dof] = typsec_parameters(typsec,rho);

    % Airfoil aerodynamic parameters 
    airfoil = 'NACA0012';
    [c_n_alpha,S1,S2,alpha1_0,alpha_ss,delta_alpha1,alpha_ds0,Df,TvL,c_n1,K0,K1,K2,x_ac,kappa,Tp,Ta,Tf0,Tv0,q0,E0,c_d0,c_m0,eta] = airfoil_parameters(airfoil,M,U,b);

    % Indicial response parameters
    [A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I] = read_indicial_params(M,beta,b,a_inf);

    % State space matrices
    [A,B] = get_SS_matrices_AEBL(U,b,beta,M,a_inf);

    %% Setup NR algorithm
    % NR convergence variables
    tol = 1e-10;        % NR tolerance 
    NR_it = 0;
    NR_max_it = 60;
    
    %% Setup Shooting method
    N = size(x0,2);
    phase_condition = 2;        % Choose phase condition - (1) for fixed state (Poincaré section), (2) for orthogonality condition
    diffs = 'f';                % Choose finite-differences scheme for df_dx calculation: 'c' for centered, 'f' for forward
    % Initialize increment vector
    y = [x0, T0]'; 
    delta_y = 10*tol; 
    eps = sqrt(delta_y'*delta_y);
    if phase_condition == 1
        % Define marker state for Poincaré section - phase of periodic solution
        marker_state = 16;
        y(marker_state) = [];
    end

    %% Run Shooting method
    while eps>tol
        [t,x] = AEBL_rkf_Sheng_Shooting([0 T0],x0,xdot0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
        f = (x(end,:)-x0)'; 
        [df_dx,df_dT] = df_dx_and_df_dT(N,T0,x0,xdot0,f,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options,diffs);
        if phase_condition == 1
            df_dx0 = df_dx;
            df_dx0(:,marker_state) = []; 
            J = [df_dx0 df_dT];    
            delta_y = -J\f;
            y = y + delta_y;
            eps = sqrt(delta_y'*delta_y);
            x0(1:marker_state-1) = y(1:marker_state-1); x0(marker_state+1:end) = y(marker_state:end-1);
            T0 = y(end);
        elseif phase_condition == 2
            f = [f; 0];
            J = [df_dx df_dT; df_dT', 0];    
            delta_y = -J\f;
            y = y + delta_y;
            eps = sqrt(delta_y'*delta_y);
            x0 = y(1:end-1)';
            T0 = y(end);
        end
        Psi = df_dx+eye(N);
        NR_it = NR_it+1;
        disp(['NR iteration ' num2str(NR_it) ': eps = ' num2str(eps)])
        [~,floquet_mtp] = eigenshuffle(Psi);
        if NR_it == NR_max_it
            warning(['Maximum number of Newton-Raphson iterations reached - stopping (eps = ' num2str(eps) ')']);
            break
        end
    end

        %% Plot Floquet multipliers evolution
        if ~exist('figure1','var')
            figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
            axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman');
            x_circle = -1:1e-3:1;
            y_circle_p = sqrt(abs(1-x_circle.^2));
            y_circle_m = -y_circle_p;
            c = jet(length(M_vec));
        end

        hold(axes1,'on');
        plot(x_circle,y_circle_p,'k-',x_circle,y_circle_m,'k-','Parent',axes1);
        hold(axes1,'on');
        plot(real(floquet_mtp),imag(floquet_mtp),'.','MarkerSize',10,'MarkerFaceColor',c(i,:),'MarkerEdgeColor',c(i,:),'LineStyle','none','Parent',axes1);
        xlabel(axes1,'$\Re(\lambda_M)$','FontWeight','normal','FontSize',axes_size);
        ylabel(axes1,'$\Im(\lambda_M)$','FontWeight','normal','FontSize',axes_size);
        title(axes1,['M = ', num2str(M)],'FontWeight','normal','FontSize',axes_size);
        grid on
        axis equal; xlim(axes1,[-1.05 1.05]); ylim(axes1,[-1.05 1.05]);
        drawnow
        hold(axes1,'off');
        
        %% Save data
        save(['typsec_Fung_M_' num2str(M,'%.6f') '_Shooting_data_new.mat'],'M','x0','xdot0','inputs_outputs_save','T0','floquet_mtp','eps','tol','NR_it','dt','options','diffs','phase_condition');
    

end