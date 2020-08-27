clc
clear all
close all
load_data = 'inicon';
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
axes_size = 20;
lw = 0.75;
ms = 5;
global divergent_flag; divergent_flag = 0;

%% Test conditions and flow
dt = 1e-5;           % Time step
tf = 300;
t_renorm = 5;
tspan = [0 tf];
h = 0e3;
[a_inf,rho,~,~] = std_atmosphere_calc(h); a_inf = 343; rho = 1.225;
M = 0.09976; 
U = M*a_inf;    
qinf = 1/2*rho*U^2;
beta = sqrt(1-M^2);

%% Choose if modifications will be allowed - 'on' for Sheng's modification, 'off' for standard BL
mods = 'on'; 

%% Airfoil structural parameters
typsec = 'mine';
[b,ah,m,I_alpha,S_alpha,Kh,Ka,Ch,Ca,Kh3,Ka3,h_F,alpha_F,plunge_dof] = typsec_parameters(typsec,rho);

%% Airfoil aerodynamic parameters 
airfoil = 'NACA0012';
[c_n_alpha,S1,S2,alpha1_0,alpha_ss,delta_alpha1,alpha_ds0,Df,TvL,c_n1,K0,K1,K2,x_ac,kappa,Tp,Ta,Tf0,Tv0,q0,E0,c_d0,c_m0,eta] = airfoil_parameters(airfoil,M,U,b);

%% Indicial response parameters
% Circulatory and non-circulatory unit step AoA indicial response - BL model
[A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I] = read_indicial_params(M,beta,b,a_inf);

%% State space matrices
[A,B] = get_SS_matrices_AEBL(U,b,beta,M,a_inf);

%% ODE solver
RKF_relax_fac_min = 1/2^0; xdotRelTol = 1e-12; xdotAbsTol = 1e-12; xdot_max_it = 10; d0 = 1e-9; r = t_renorm/dt;
options = {'hlim',dt,'relax_fac_min',RKF_relax_fac_min,'xdotRelTol',xdotRelTol,'xdotAbsTol',xdotAbsTol,'xdot_max_it',xdot_max_it,'d0',d0,'r',r};

%% Initial conditions 
typsec_filename = get_typsec_filename(typsec,b,ah,m,I_alpha,S_alpha,Kh,Ka,Ka3,h_F,alpha_F,rho,h);
if strcmp(load_data,'inicon')
    filename = [typsec_filename,'_M_', num2str(M,'%.5f'), '.mat'];
    load(filename,'x0_save','xdot0_save','inputs_outputs_save')
    x0 = x0_save;
    xdot0 = xdot0_save;
    inputs_outputs = inputs_outputs_save;
    renorms_total = [];
    L1_total = [];
    norm_d_total = [];
elseif strcmp(load_data,'continue')
    filename = [typsec_filename,'_M_', num2str(M,'%.5f') '_dt_' num2str(dt,'%.0e') '_tr_' num2str(t_renorm) '_Lyapunov_2PM.mat']; 
    load(filename,'x0_save','tf_save','inputs_outputs_save','L1')
    x0 = x0_save;
    xdot0 = zeros(1,length(x0));
    inputs_outputs = inputs_outputs_save;
    tspan = [tf_save tf_save+tf];
end

%% Calculate flow
if strcmp(mods,'off') % Standard BL
%     [t,x,y,xdot] = AEBL_rkf_St(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,mods,q0,options);
else                  % Sheng's modifications
    [t,d,x0_save,inputs_outputs_save] = AEBL_rkf_Sheng_Lyapunov_2PM(tspan,x0,xdot0,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
end

%% Calculate disturbance evolution and mLE    
% Maximal Lyapunov exponent
n = length(t);
if r < n
    p = round(n/r);
    L1 = zeros(p,1);
    norm_d = zeros(p,1);
    cum_sum_log = zeros(p,1);
    for i=1:p
        norm_d(i) = norm(d(i*round(r),:));
        if i == 1
            cum_sum_log(i) = log(norm_d(i)/d0);
        else
            cum_sum_log(i) = cum_sum_log(i-1) + log(norm_d(i)/d0);
        end
        L1(i) = 1/t(i*round(r))*cum_sum_log(i);
    end
else
    p = 0;
    norm_d = d0*ones(n,1);
    L1 = zeros(n,1);
    for i=1:n-1
        norm_d(i+1) = norm(d(i+1,:));
        L1(i+1) = 1/t(i+1)*log(norm_d(i+1)/d0);
    end
end
mLE = L1(end) % maximal Lyapunov Exponent

%% Increment total vectors
if ~isempty(renorms_total)
    r1 = 1;
else
    r1 = length(renorms_total);
end
renorms = r1:p+r1-1;
renorms_total = [renorms_total renorms];
L1_total = [L1_total; L1];
norm_d_total = [norm_d_total; norm_d];

%% Plots
if r < n
    figure1 = figure('InvertHardcopy','off','Color',[1 1 1],'Units','normalized','Position',[0 0.28 1 0.58]);
    axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman','XScale','linear','YScale','log');
    hold(axes1,'on'); 
    semilogy(renorms_total,norm_d_total/d0,'k-o','LineWidth',lw,'MarkerSize',ms,'Parent',axes1);
    xlabel('Renormalizations','FontWeight','normal','FontSize',axes_size);
    ylabel('$\|d\|/\|d_0\|$','FontWeight','normal','FontSize',axes_size);
    title(['M = ', num2str(M), ', dt = ', num2str(dt,'%.0e'), ', t = ' num2str(t(end))],'FontWeight','normal','FontSize',axes_size);
    grid on

    figure2 = figure('InvertHardcopy','off','Color',[1 1 1],'Units','normalized','Position',[0 0.28 1 0.58]);
    axes2 = axes('Parent',figure2,'FontSize',axes_size,'FontName','times new roman','XScale','log','YScale','log');
    hold(axes2,'on'); 
    loglog(renorms_total,L1_total,'k-o','LineWidth',lw,'Parent',axes2);
    xlabel('Renormalizations','FontWeight','normal','FontSize',axes_size);
    ylabel('$\hat{L}_1$','FontWeight','normal','FontSize',axes_size);
    title(['M = ', num2str(M), ', dt = ', num2str(dt,'%.0e'), ', t = ' num2str(t(end))],'FontWeight','normal','FontSize',axes_size);
    grid on
else
    figure1 = figure('InvertHardcopy','off','Color',[1 1 1],'Units','normalized','Position',[0 0.28 1 0.58]);
    axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman','XScale','linear','YScale','log');
    hold(axes1,'on'); 
    semilogy(t,norm_d/d0,'k-','LineWidth',lw,'Parent',axes1);
    xlabel('t [s]','FontWeight','normal','FontSize',axes_size);
    ylabel('$\|d\|/\|d_0\|$','FontWeight','normal','FontSize',axes_size);
    title(['M = ', num2str(M), ', dt = ', num2str(dt,'%.0e'), ', t = ' num2str(t(end))],'FontWeight','normal','FontSize',axes_size);
    grid on

    figure2 = figure('InvertHardcopy','off','Color',[1 1 1],'Units','normalized','Position',[0 0.28 1 0.58]);
    axes2 = axes('Parent',figure2,'FontSize',axes_size,'FontName','times new roman','XScale','log','YScale','log');
    hold(axes2,'on'); 
    loglog(t,L1,'k-','LineWidth',lw,'Parent',axes2);
    xlabel('t [s]','FontWeight','normal','FontSize',axes_size);
    ylabel('$\hat{L}_1$','FontWeight','normal','FontSize',axes_size);
    title(['M = ', num2str(M), ', dt = ', num2str(dt,'%.0e'), ', t = ' num2str(t(end))],'FontWeight','normal','FontSize',axes_size);
    grid on
end


%% Save data
tf_save = t(end);
if r == 1 
    filename = [typsec_filename,'_M_', num2str(M,'%.5f') '_Lyapunov_2PM_dt_' num2str(dt,'%.0e') '_tr_dt.mat'];
elseif r < n
    filename = [typsec_filename,'_M_', num2str(M,'%.5f') '_Lyapunov_2PM_dt_' num2str(dt,'%.0e') '_tr_' num2str(t_renorm) '.mat']; 
else
    filename = [typsec_filename,'_M_', num2str(M,'%.5f') '_Lyapunov_2PM_dt_' num2str(dt,'%.0e') '_tr_inf.mat']; 
end
save(filename,'mLE','dt','d0','renorms_total','L1_total','norm_d_total','x0_save','tf_save','inputs_outputs_save')
