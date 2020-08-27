clc
clear all
close all
set(0,'DefaultTextInterpreter','latex')
axes_size = 20;
lw = 0.75;

global jM continue_run inputs_outputs_save; jM = 1; continue_run = 'n'; 
M_init = 0.07651;
M_final = 0.07651;
dM = 1e-4;
M_vec = M_init:dM:M_final;
t_discard = 20;

for i=1:length(M_vec)
M = M_vec(i);
%% Test conditions and flow
dt = 5e-5;
tf = 120; % Runtime 
h = 0; [a_inf,rho,~,~] = std_atmosphere_calc(h); a_inf = 343; rho = 1.225;
U = M*a_inf;    
qinf = 1/2*rho*U^2;
beta = sqrt(1-M^2);

%% Choose if modifications will be allowed - 'on' for Sheng's modification, 'off' for standard BL
mods = 'on'; 

%% Airfoil structural parameters
typsec = 'Fung';
[b,ah,m,I_alpha,S_alpha,Kh,Ka,Ch,Ca,Kh3,Ka3,h_F,alpha_F,plunge_dof] = typsec_parameters(typsec,rho);
typsec_filename = get_typsec_filename(typsec,b,ah,m,I_alpha,S_alpha,Kh,Ka,Ka3,h_F,alpha_F,rho,h);

%% Airfoil aerodynamic parameters 
airfoil = 'NACA0012';
[c_n_alpha,S1,S2,alpha1_0,alpha_ss,delta_alpha1,alpha_ds0,Df,TvL,c_n1,K0,K1,K2,x_ac,kappa,Tp,Ta,Tf0,Tv0,q0,E0,c_d0,c_m0,eta] = airfoil_parameters(airfoil,M,U,b);

%% Indicial response parameters
% Circulatory and non-circulatory unit step AoA indicial response - BL model
[A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I] = read_indicial_params(M,beta,b,a_inf);

%% State space matrices
[A,B] = get_SS_matrices_AEBL(U,b,beta,M,a_inf);

%% Load initial conditions 
load_from_file = 'y';
if strcmp(load_from_file,'y')
    filename = ['typsec_Fung_M_' num2str(M,'%.5f') '.mat'];
    load(filename,'tf_save','x0_save','xdot0_save','y0_save','inputs_outputs_save');
    ti = tf_save; 
    tf = ti+tf;
    x0 = x0_save;
    y0 = y0_save;
    xdot0 = xdot0_save;
    tspan = [ti tf];
else
    h_0 = 0.05*b; 
    h_dot0 = 0; 
    h_2dot0 = 0; 
    alpha_0 = 25*pi/180; 
    alpha_dot0 = 0; 
    alpha_2dot0 = 0;
    % Initialize BL and structural states 
    [x0,y0,xdot0] = initialize_BL(h_0,h_dot0,h_2dot0,alpha_0,alpha_dot0,alpha_2dot0,A,B,U,M,b,ah,alpha1_0,delta_alpha1,S1,S2,c_n1,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,eta,E0,Df,alpha_ds0,alpha_ss,q0,mods);
    ti = 0;
    tspan = [ti tf];
end

%% ODE solver
RKF_relax_fac_min = 1/2^0; RKFtol = 1e-10; xdotRelTol = 1e-12; xdotAbsTol = 1e-12; xdot_max_it = 10; disp_boundaries = 1; 
options = {'SA_LCO_detect',0,'MA_LCO_detect',0,'hlim',dt,'relax_fac_min',RKF_relax_fac_min,'RKFtol',RKFtol,'xdotRelTol',xdotRelTol,'xdotAbsTol',xdotAbsTol,'xdot_max_it',xdot_max_it,'disp_boundaries',disp_boundaries};
if strcmp(mods,'off') % Standard BL
    [t, x, y, xdot] = AEBL_rkf_St(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,mods,q0,options);
else                  % Sheng's modifications
    [t, x, y, xdot] = AEBL_rkf_Sheng(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
end

%% Output variables 
xi = x(:,15)/b;
alpha = x(:,16); 
xi_dot = x(:,13)/b;
alpha_dot = x(:,14);
alpha_bar = y(:,4);
alpha_bar_dot = y(:,5);
alpha_cr = y(:,6);
x9 = x(:,9);
x10 = x(:,10);

%% Poincaré sections
% Discard possible transients
t_red = t-t(1);
p = find(t_red<=t_discard,1,'last');
alpha_red = alpha(p:end); alpha_dot_red = alpha_dot(p:end); xi_red = xi(p:end); xi_dot_red = xi_dot(p:end); alpha_bar_red = alpha_bar(p:end); alpha_bar_dot_red = alpha_bar_dot(p:end);
[alpha_P,alpha_dot_P,xi_P,xi_dot_P,alpha_bar_P,alpha_bar_dot_P] = poincare_section(alpha_red,alpha_dot_red,xi_red,xi_dot_red,alpha_bar_red,alpha_bar_dot_red);

%% Outputs spectra
[f_vec,SP1_alpha,SP1_xi,SP1_alpha_bar] = get_FFT(t,alpha,xi,alpha_bar,dt);

%% Save variables 
filename = [typsec_filename,'_M_', num2str(M,'%.5f'), '_spectrum_and_poincare.mat'];
save (filename,'M','dt','f_vec','t','xi','xi_dot','alpha','alpha_dot','alpha_bar','alpha_bar_dot','alpha_cr','x9','x10','SP1_alpha','SP1_xi','SP1_alpha_bar','alpha_P','alpha_dot_P','xi_P','xi_dot_P','alpha_bar_P','alpha_bar_dot_P'); 

%% Plots
close all

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman');
hold(axes1,'on'); 
plot(alpha_P*180/pi,alpha_dot_P*180/pi,'k.','LineWidth',lw,'Parent',axes1);
xlabel('$\alpha$ [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('$\dot{\alpha}$ [deg/s]','FontWeight','normal','FontSize',axes_size);
title(['M = ', num2str(M)],'FontWeight','normal','FontSize',axes_size);
grid off
drawnow

figure2 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes2 = axes('Parent',figure2,'FontSize',axes_size,'FontName','times new roman');
hold(axes2,'on'); 
plot(xi_P,xi_dot_P,'k.','LineWidth',lw,'Parent',axes2);
xlabel('$\xi$','FontWeight','normal','FontSize',axes_size);
ylabel('$\dot{\xi}$ [/s]','FontWeight','normal','FontSize',axes_size);
title(['M = ', num2str(M)],'FontWeight','normal','FontSize',axes_size);
grid off
drawnow

figure3 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes3 = axes('Parent',figure3,'FontSize',axes_size,'FontName','times new roman','YScale','log','XScale','log');
hold(axes3,'on');
loglog(f_vec*2*pi*b/U,SP1_xi,'k-','LineWidth',lw,'Parent',axes3);
xlabel('$k$','FontWeight','normal','FontSize',axes_size);
ylabel('Amplitude spectrum of $\xi$','FontWeight','normal','FontSize',axes_size);
title(['M = ', num2str(M)],'FontWeight','normal','FontSize',axes_size);
grid off
drawnow

figure4 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes4 = axes('Parent',figure4,'FontSize',axes_size,'FontName','times new roman','YScale','log','XScale','log');
hold(axes4,'on');
loglog(f_vec*2*pi*b/U,SP1_alpha,'k-','LineWidth',lw,'Parent',axes4);
xlabel('$k$','FontWeight','normal','FontSize',axes_size);
ylabel('Amplitude spectrum of $\alpha$','FontWeight','normal','FontSize',axes_size);
title(['M = ', num2str(M)],'FontWeight','normal','FontSize',axes_size);
grid off
drawnow

end