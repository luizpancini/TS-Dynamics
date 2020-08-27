clc
clearvars -except xdot0_save x0_save tf_save
close all

global continue_run jM; jM = 0;
i = 1;
i_final = 10; % Number of times to run (multiply by tf to get total simulation time)
alpha_P_total = [];
xi_P_total = [];
alpha_dot_P_total = [];
xi_dot_P_total = [];

while i<=i_final
clearvars -except i i_final continue_run jM xdot0_save x0_save tf_save alpha_P_total xi_P_total alpha_dot_P_total xi_dot_P_total 
jM = jM+1;
i
%% Test conditions and flow
tf = 1; % maximum single runtime
if jM == 1, ti = 0; else, ti = tf_save; end
dt = 2e-4;
tspan = [ti tf+ti];
h = 0; [a_inf,rho,~,~] = std_atmosphere_calc(h);
M = 0.0770; 
U = M*a_inf;    
qinf = 1/2*rho*U^2;
beta = sqrt(1-M^2);

%% Choose if modifications will be allowed - 'on' for Sheng's modification, 'off' for standard BL
mods = 'on'; 

%% Airfoil structural parameters
typsec = 'Fung';
[b,ah,m,I_alpha,S_alpha,Kh,Ka,Ch,Ca,Kh3,Ka3,h_F,alpha_F,plunge_dof] = typsec_parameters(typsec,rho);

%% Airfoil aerodynamic parameters 
airfoil = 'NACA0012';
[c_n_alpha,S1,S2,alpha1_0,alpha_ss,delta_alpha1,alpha_ds0,Df,TvL,c_n1,K0,K1,K2,x_ac,kappa,Tp,Ta,Tf0,Tv0,q0,E0,c_d0,c_m0,eta] = airfoil_parameters(airfoil,M,U,b);

%% Initial conditions 
h_0 = 0*b; 
h_dot0 = 0; 
h_2dot0 = 0; 
alpha_0 = 20*pi/180; 
alpha_dot0 = 0; 
alpha_2dot0 = 0;

%% Indicial response parameters
% Circulatory and non-circulatory unit step AoA indicial response - BL model
[A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I] = read_indicial_params(M,beta,b,a_inf);

%% State space matrices
[A,B,~,~,A_bar,SM_inv] = get_SS_matrices_AEBLgust(U,b,beta,M,a_inf,c_n_alpha,m,S_alpha,I_alpha,Ch,Ca,Kh,Ka);

%% Initialize BL and structural states 
[x0,y0,xdot0] = initialize_BL(h_0,h_dot0,h_2dot0,alpha_0,alpha_dot0,alpha_2dot0,A,B,U,M,b,ah,alpha1_0,delta_alpha1,S1,S2,c_n1,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,eta,E0,Df,alpha_ds0,alpha_ss,q0,mods);

%% Initialize from previous run
continue_run = 'y'; % 'y' for yes and 'n' for no
if strcmp(continue_run,'y') && jM>1
    x0 = x0_save;
    xdot0 = xdot0_save;
end

%% ODE solver
if strcmp(mods,'off') % Standard BL
    [t, x, ~, xdot] = AEBL_rkf_St(dt,tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,mods,q0);
else                  % Sheng's modifications
    [t, x, ~, xdot] = AEBL_rkf_Sheng(dt,tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1);
end

% Truncate pre-allocated vectors
global final_it;
t = t(1:final_it); x = x(1:final_it,:); xdot = xdot(1:final_it,:);

%% Output variables 
xi_dot = x(:,13)/b;
alpha_dot = x(:,14);
xi = x(:,15)/b;
alpha = x(:,16);

[alpha_P,alpha_dot_P,xi_P,xi_dot_P] = poincare_section(alpha,alpha_dot,xi,xi_dot);

%% Save to total vectors
alpha_P_total = [alpha_P_total alpha_P];
xi_P_total = [xi_P_total xi_P];
alpha_dot_P_total = [alpha_dot_P_total alpha_dot_P];
xi_dot_P_total = [xi_dot_P_total xi_dot_P];

%% Save final conditions for next run's initial conditions 
x0_save = x(end,:);
xdot0_save = xdot(end,:);
tf_save = tspan(end);

i = i+1;
end

%% Save total vectors
typsec_filename = get_typsec_filename(typsec,b,ah,m,I_alpha,S_alpha,Kh,Ka,Ka3,h_F,alpha_F,rho);
filename = [typsec_filename,'_M_', num2str(M,'%.5f'),'_poincare.mat'];
save (filename,'M','alpha_P_total','xi_P_total','alpha_dot_P_total','xi_dot_P_total'); 

%% Final plots
set(0,'DefaultTextInterpreter','latex')
axes_size = 20;
lw = 0.75;

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman');
hold(axes1,'on'); 
plot(alpha_P_total*180/pi,alpha_dot_P_total*180/pi,'k.','Parent',axes1);
xlabel('$\alpha$ [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('$\dot{\alpha}$ [deg/s]','FontWeight','normal','FontSize',axes_size);
grid off

figure2 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes2 = axes('Parent',figure2,'FontSize',axes_size,'FontName','times new roman');
hold(axes2,'on'); 
plot(xi_P_total,xi_dot_P_total,'k.','Parent',axes2);
xlabel('$\xi$','FontWeight','normal','FontSize',axes_size);
ylabel('$\dot{\xi}$ [/s]','FontWeight','normal','FontSize',axes_size);
grid off
