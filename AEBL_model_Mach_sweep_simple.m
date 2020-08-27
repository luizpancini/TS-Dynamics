clc
clear all
close all

%% TS known data
alpha_I = 0*pi/180; % Wind-off angle (fixed point)

%% Mach sweep options
M = 0.07693;
dM = 1e-5;

%% Run the Mach sweep
fixed_point_flag = 0;
global divergent_flag continue_run jM; divergent_flag = 0; continue_run = 'y'; jM = 0;

while fixed_point_flag == 0

%% Test conditions and flow
dt = 5e-5;
tf = 60; % Runtime 
h = 0e3; [a_inf,rho,~,~] = std_atmosphere_calc(h); a_inf = 343; rho = 1.225;
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

%% Initial conditions
load(['typsec_Fung_M_',num2str(M,'%.6f'),'_LCO_data.mat'], 'x0_save','xdot0_save','inputs_outputs_save','tf_save');
x0 = x0_save;
xdot0 = xdot0_save;
y0 = zeros(1,20);
tspan = [tf_save tf_save+tf];

%% Indicial response parameters
% Circulatory and non-circulatory unit step AoA indicial response - BL model
[A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I] = read_indicial_params(M,beta,b,a_inf);

%% State space matrices
[A,B] = get_SS_matrices_AEBL(U,b,beta,M,a_inf);

%% Initialize BL and structural states - only y0
% [~,y0,~] = initialize_BL(h_0,h_dot0,h_2dot0,alpha_0,alpha_dot0,alpha_2dot0,A,B,U,M,b,ah,alpha1_0,delta_alpha1,S1,S2,c_n1,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,eta,E0,Df,alpha_ds0,alpha_ss,q0,mods);

%% ODE solver
options = {'hlim',dt,'SA_LCO_detect',1,'MA_LCO_detect',0,'relax_fac_min',1,'xdotRelTol',1e-12,'xdotAbsTol',1e-12,'RKFtol',1e-10};
if strcmp(mods,'off') % Standard BL
    [t, x, y, xdot] = AEBL_rkf_St(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,mods,q0,options);
else                  % Sheng's modifications
    [t, x, y, xdot] = AEBL_rkf_Sheng(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
end

xi = x(:,15)/b;
alpha = x(:,16);

%% Plots
set(0,'DefaultTextInterpreter','latex')
axes_size = 20;
lw = 0.75;
close all

%% Structure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman');
hold(axes1,'on'); 
plot(t,alpha*180/pi,'k-','LineWidth',lw,'Parent',axes1);
xlabel('t [s]','FontWeight','normal','FontSize',axes_size);
ylabel('$\alpha$ [deg]','FontWeight','normal','FontSize',axes_size);
grid off
drawnow

figure2 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes2 = axes('Parent',figure2,'FontSize',axes_size,'FontName','times new roman');
hold(axes2,'on'); 
plot(t,xi,'k-','LineWidth',lw,'Parent',axes2);
xlabel('t [s]','FontWeight','normal','FontSize',axes_size);
ylabel('$\xi$','FontWeight','normal','FontSize',axes_size);
grid off
drawnow

%% PSD and FFT 
[f_fft,SP1_alpha,SP1_xi,alpha_PSD,xi_PSD,~,~,N0_alpha,N0_xi] = get_PSD_and_FFT(t,alpha,xi,dt);
% 
% figure16 = figure('InvertHardcopy','off','Color',[1 1 1]);
% axes16 = axes('Parent',figure16,'FontSize',axes_size,'FontName','times new roman');
% loglog(f_fft,alpha_PSD,'k-','LineWidth',lw,'Parent',axes16);
% xlabel('f [Hz]','FontWeight','normal','FontSize',axes_size);
% ylabel('$\Phi(\alpha)$ [rad$^2$/Hz]','FontWeight','normal','FontSize',axes_size);
% grid off
% 
% figure17 = figure('InvertHardcopy','off','Color',[1 1 1]);
% axes17 = axes('Parent',figure17,'FontSize',axes_size,'FontName','times new roman');
% loglog(f_fft,xi_PSD,'k-','LineWidth',lw,'Parent',axes17);
% xlabel('f [Hz]','FontWeight','normal','FontSize',axes_size);
% ylabel('$\Phi(\xi)$ [rad$^2$/Hz]','FontWeight','normal','FontSize',axes_size);
% grid off

figure25 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes25 = axes('Parent',figure25,'FontSize',axes_size,'FontName','times new roman','YScale','log','XScale','log');
loglog(f_fft,SP1_xi,'k-','LineWidth',lw,'Parent',axes25);
xlabel('f [Hz]','FontWeight','normal','FontSize',axes_size);
ylabel('Amplitude spectrum of $\xi$','FontWeight','normal','FontSize',axes_size);
grid off
drawnow

figure26 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes26 = axes('Parent',figure26,'FontSize',axes_size,'FontName','times new roman','YScale','log','XScale','log');
loglog(f_fft,SP1_alpha,'k-','LineWidth',lw,'Parent',axes26);
xlabel('f [Hz]','FontWeight','normal','FontSize',axes_size);
ylabel('Amplitude spectrum of $\alpha$','FontWeight','normal','FontSize',axes_size);
grid off
drawnow

%% Save final conditions for next run's initial conditions
x0_save = x(end,:);
y0_save = y(end,:);
xdot0_save = xdot(end,:);
tf_save = t(end);

%% Save data
global t_alpha_lpeak t_alpha_upeak inputs_outputs_save
LCO_period = mean([diff(t_alpha_lpeak), diff(t_alpha_upeak)]);
save(['typsec_Fung_M_',num2str(M,'%.6f'),'_LCO_data.mat'], 'x0_save','xdot0_save','inputs_outputs_save','LCO_period','tf_save');
 
%% Check for stationary solution
alpha_max = max(alpha)*180/pi
if abs(alpha_I-alpha_max) < 1e-2
    warning('Oscillations too close to wind-off angle reached - assuming fixed point has been reached and stopping')
    fixed_point_flag = 1;
    break;
end

M = M+dM;
end
