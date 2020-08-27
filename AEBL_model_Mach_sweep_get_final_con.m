clc
clear all
close all

%% TS known data
M_F = 0.07693;      % Flutter Mach speed
alpha_I = 0*pi/180; % Wind-off angle (fixed point)

%% Mach sweep options
Mach_init = 0.07600;
dM = -1e-5;

%% Run the Mach sweep
M = Mach_init; LCO_flag = 0; fixed_point_flag = 0;
global jM divergent_flag continue_run LCO_reached inputs_outputs_save; jM = 0; divergent_flag = 0;

while fixed_point_flag == 0 && divergent_flag == 0
jM = jM+1;
if LCO_flag == 1 % if LCO has been reached, advance the Mach
    M = M+dM;
    LCO_reached = [];
end
M

%% Test conditions and flow
dt = 2e-4;
tf = 10; % Maximum single runtime 
if jM == 1, ti = 0; else, ti = tf_save; end
tspan = [ti tf];
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

%% Initial conditions 
h_0 = 0.1*b; 
h_dot0 = 0; 
h_2dot0 = 0; 
alpha_0 = 25*pi/180; 
alpha_dot0 = 0; 
alpha_2dot0 = 0;

%% Indicial response parameters
% Circulatory and non-circulatory unit step AoA indicial response - BL model
[A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I] = read_indicial_params(M,beta,b,a_inf);

%% State space matrices
[A,B] = get_SS_matrices_AEBL(U,b,beta,M,a_inf);

%% Initialize BL and structural states 
[x0,y0,xdot0] = initialize_BL(h_0,h_dot0,h_2dot0,alpha_0,alpha_dot0,alpha_2dot0,A,B,U,M,b,ah,alpha1_0,delta_alpha1,S1,S2,c_n1,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,eta,E0,Df,alpha_ds0,alpha_ss,q0,mods);

%% Initialize from previous run
continue_run = 'y'; % 'y' for yes and 'n' for no
if strcmp(continue_run,'y') && jM>1
    x0 = x0_save;
    xdot0 = xdot0_save;
    tf = tf+tf_save;
    tspan = [tf_save tf];
end

%% ODE solver
options = {'hlim',dt,'SA_LCO_detect',1,'MA_LCO_detect',1,'relax_fac_min',1};
if strcmp(mods,'off') % Standard BL
    [t, x, y, xdot] = AEBL_rkf_St(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,delta_alpha1,S1,S2,c_n1,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Tp,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,mods,q0,options);
else                  % Sheng's modifications
    [t, x, y, xdot] = AEBL_rkf_Sheng(tspan,x0,xdot0,y0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
end

alpha = x(:,16);

%% Plots
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
axes_size = 20;
lw = 0.75;
close all   
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman');
hold(axes1,'on'); 
plot(t,alpha*180/pi,'k-','LineWidth',lw,'Parent',axes1);
ylabel('$\alpha$ [deg]','FontWeight','normal','FontSize',axes_size);
xlabel('t [s]','FontWeight','normal','FontSize',axes_size);
grid on
title(['M = ', num2str(M)],'FontWeight','normal','FontSize',axes_size);
drawnow 

%% Check if fixed point has been reached
alpha_max = max(alpha)*180/pi;
if abs(alpha_I-alpha_max) < 1e-2
    warning('Oscillations too close to wind-off angle reached - assuming fixed point has been reached and stopping')
    fixed_point_flag = 1;
    break;
end

%% Save final conditions for next run's initial conditions
x0_save = x(end,:);
y0_save = y(end,:);
xdot0_save = xdot(end,:);
tf_save = tf;

% Save 
if ~isempty(LCO_reached) && divergent_flag == 0
    LCO_flag = 1; 
    inputs_outputs0_save = inputs_outputs_save;
    filename = [typsec_filename,'_M_', num2str(M,'%.5f'), '_final_con.mat'];
    typsec_filepath = filename;
    save (typsec_filepath,'M_F','M','x0_save','y0_save','xdot0_save','inputs_outputs0_save','a_inf','rho')
end 

end