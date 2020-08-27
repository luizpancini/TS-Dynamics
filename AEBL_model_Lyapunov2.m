clc
clear all
close all
load_inicon_data = 'y';
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
axes_size = 20;
lw = 0.75;
ms = 5;
global divergent_flag; divergent_flag = 0;

%% Test conditions and flow
dt = 1e-5;           % Time step
tf = 300;
t_renorm = 1e9;
tspan = [0 tf];
h = 0e3;
[a_inf,rho,~,~] = std_atmosphere_calc(h); a_inf = 343; rho = 1.225;
M = 0.09577; 
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
 
%% Initial conditions 
if strcmp(load_inicon_data,'y') 
    typsec_filename = get_typsec_filename(typsec,b,ah,m,I_alpha,S_alpha,Kh,Ka,Ka3,h_F,alpha_F,rho,h);
    filename = [typsec_filename,'_M_', num2str(M,'%.5f'), '.mat']; 
    load(filename,'x0_save','xdot0_save','inputs_outputs_save')
    x0 = x0_save;
    xdot0 = xdot0_save;
    inputs_outputs = inputs_outputs_save;
else
    h_0 = 0.05*b; 
    h_dot0 = 0; 
    h_2dot0 = 0; 
    alpha_0 = 25*pi/180; 
    alpha_dot0 = 0; 
    alpha_2dot0 = 0;
    % Initialize BL and structural states 
    [x0,y0,xdot0] = initialize_BL(h_0,h_dot0,h_2dot0,alpha_0,alpha_dot0,alpha_2dot0,A,B,U,M,b,ah,alpha1_0,delta_alpha1,S1,S2,c_n1,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,eta,E0,Df,alpha_ds0,alpha_ss,q0,mods);
    inputs_outputs = zeros(7,2);
end

%% Set disturbance 
d0 = 1;
d_i = d0*ones(length(x0),1); d_i = d0*d_i/norm(d_i); 

%% Calculate flow
RKF_relax_fac_min = 1/2^0; xdotRelTol = 1e-12; xdotAbsTol = 1e-12; xdot_max_it = 10; FD_method = 'c2'; int_method = 1; r = t_renorm/dt;
options = {'hlim',dt,'relax_fac_min',RKF_relax_fac_min,'xdotRelTol',xdotRelTol,'xdotAbsTol',xdotAbsTol,'xdot_max_it',xdot_max_it,'FD_method',FD_method,'d0',d0,'r',r,'int_method',int_method};
if strcmp(mods,'off') % Standard BL
%     [t,x,y,xdot] = AEBL_rkf_St(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,mods,q0,options);
else                  % Sheng's modifications
    [t,d,x0_save,inputs_outputs] = AEBL_rkf_Sheng_Lyapunov(tspan,x0,xdot0,inputs_outputs,d_i,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
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
        norm_d(i) = norm(d(:,i*round(r)));
        if i == 1
            cum_sum_log(i) = log(norm_d(i)/d0);
        else
            cum_sum_log(i) = cum_sum_log(i-1) + log(norm_d(i)/d0);
        end
        L1(i) = 1/t(i*round(r))*cum_sum_log(i);
    end
else
    norm_d = d0*ones(n,1);
    L1 = zeros(n,1);
    for i=1:n-1
        norm_d(i+1) = norm(d(:,i+1));
        L1(i+1) = 1/t(i+1)*log(norm_d(i+1)/d0);
    end
end

mLE = L1(end)
        
%% Plots
if r < n
    figure1 = figure('InvertHardcopy','off','Color',[1 1 1],'Units','normalized','Position',[0 0.28 1 0.58]);
    axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman','XScale','log','YScale','log');
    hold(axes1,'on'); 
    loglog(renorms,L1,'k-o','LineWidth',lw,'Parent',axes1);
    xlabel('Renormalizations','FontWeight','normal','FontSize',axes_size);
    ylabel('$\hat{L}_1$','FontWeight','normal','FontSize',axes_size);
    title(['M = ', num2str(M), ', dt = ', num2str(dt,'%.0e'), ', t = ' num2str(t(end))],'FontWeight','normal','FontSize',axes_size);
    grid on
else
    figure1 = figure('InvertHardcopy','off','Color',[1 1 1],'Units','normalized','Position',[0 0.28 1 0.58]);
    axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman','XScale','log','YScale','log');
    hold(axes1,'on'); 
    loglog(t,abs(L1),'k-','LineWidth',lw,'Parent',axes1);
    xlabel('t [s]','FontWeight','normal','FontSize',axes_size);
    ylabel('$\hat{L}_1$','FontWeight','normal','FontSize',axes_size);
    title(['M = ', num2str(M), ', dt = ', num2str(dt,'%.0e'), ', t = ' num2str(t(end))],'FontWeight','normal','FontSize',axes_size);
    grid on
end

%% Save final conditions for next run's initial conditions
d_save = d(:,end);
tf_save = tspan(end);
inputs_outputs_save = inputs_outputs;

%% Save data
typsec_filename = get_typsec_filename(typsec,b,ah,m,I_alpha,S_alpha,Kh,Ka,Ka3,h_F,alpha_F,rho,h);
if r == 1
    filename = [typsec_filename,'_M_' num2str(M,'%.5f') '_Lyapunov_dt_' num2str(dt,'%.0e') '_tr_dt.mat'];
elseif r < n
    filename = [typsec_filename,'_M_' num2str(M,'%.5f') '_Lyapunov_dt_' num2str(dt,'%.0e') '_tr_' num2str(t_renorm) '.mat'];     
else
    filename = [typsec_filename,'_M_' num2str(M,'%.5f') '_Lyapunov_dt_' num2str(dt,'%.0e') '_tr_inf.mat']; 
end
save(filename,'mLE','L1','t','norm_d','x0_save','tf_save','inputs_outputs_save','d_save','dt','FD_method','d0','t_renorm','r')
 