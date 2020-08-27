function [c_n_alpha,S1,S2,alpha1_0,alpha_ss,delta_alpha1,alpha_ds0,Df,TvL,c_n1,K0,K1,K2,x_ac,kappa,Tp,Ta,Tf0,Tv0,q0,E0,c_d0,c_m0,eta] = airfoil_parameters(airfoil,M,U,b)

%% Airfoil parameters table - as functions of Mach
if strcmp(airfoil,'NACA0012')
    Mach_range =                [0.035;   0.072;   0.110;   0.185;   0.215;   0.25;    0.3;     0.4;     0.5;     0.6;     0.7;     0.75;    0.8];
    c_n_alpha_range =    180/pi*[0.105;   0.108;   0.108;   0.110;   0.113;   0.115;   0.116;   0.113;   0.117;   0.127;   0.154;   0.175;   0.216];
    S1_range =           pi/180*[3.0;     3.0;     3.0;     3.5;     3.5;     3.5;     3.0;     3.25;    3.5;     4.0;     4.5;     3.5;     0.70];
    S2_range =           pi/180*[1.5;     1.5;     1.5;     2.0;     2.0;     2.0;     1.5;     1.6;     1.2;     0.7;     0.5;     0.8;     0.18];
    alpha1_0_range =     pi/180*[12.4;    13.8;    15.2;    16.3;    16.5;    16.0;    13.7;    12.5;    10.5;    8.5;     5.6;     3.5;     0.7];
    alpha_ss_range =     pi/180*[11.8;    13.5;    14.6;    16.0;    16.0;    16.0;    13.9;    12.5;    10.5;    8.5;     5.6;     3.5;     0.7];
    delta_alpha1_range = pi/180*[2.0;     1.0;     2.5;     2.5;     2.5;     2.5;     0.5;     2.0;     1.45;    1.0;     0.8;     0.2;     0.1];
    alpha_ds0_range =    pi/180*[13.5;    17.5;    18.5;    19.0;    17.0;    15.0;    14.0;    13.5;    13.5;    13.5;    13.5;    13.5;    13.5];
    Df_range =                  [8.0;     8.0;     8.0;     8.0;     8.0;     8.0;     8.0;     7.75;    6.2;     6.0;     5.9;     5.5;     4.0];
    TvL_range =                 [4.0;     5.0;     5.0;     5.0;     5.0;     5.0;     5.0;     9.0;     9.0;     9.0;     9.0;     9.0;     9.0];
    c_n1_range =                [1.25;    1.45;    1.60;    1.80;    1.80;    1.85;    1.60;    1.2;     1.05;    0.92;    0.68;    0.5;     0.18];
    K0_range =                  [0.0025;  0.0025;  0.0025;  0.0025;  0.0025;  0.0025;  0.0025;  0.006;   0.02;    0.038;   0.03;    0.001;   -0.01];
    K1_range =                  [-0.120;  -0.120;  -0.120;  -0.120;  -0.120;  -0.120;  -0.120;  -0.135;  -0.125;  -0.12;   -0.09;   -0.13;   0.02];
    K2_range =                  [0.04;    0.04;    0.04;    0.04;    0.04;    0.04;    0.04;    0.05;    0.04;    0.04;    0.15;    -0.02;   -0.01];
    kappa_range =               [2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0];
    Tp_range =                  [1.7;     1.7;     1.7;     1.7;     1.7;     1.7;     1.7;     1.8;     2.0;     2.5;     3.0;     3.3;     4.3];
    Ta_range =                  [3.9;     3.9;     3.9;     3.9;     3.9;     3.9;     3.9;     3.9;     3.9;     3.9;     3.9;     3.9;     3.9];
    Tf0_range =                 [3.0;     3.0;     3.0;     3.0;     3.0;     3.0;     3.0;     2.5;     2.2;     2.0;     2.0;     2.0;     2.0];
    Tv0_range =                 [6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     4.0];
    q0_range =                  [0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.010];
    c_m0_range =                [0.000;   0.000;   0.000;   0.000;   0.000;   0.000;   0.000;   0.000;   0.000;   0.000;   0.000;   0.000;   0.000];
    c_d0_range =                [0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.008;   0.0077;  0.0078;  0.0078;  0.0079;  0.0114];
    eta_range =                 [0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95];
    E0_range =                  [0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00];
else
    error('airfoil not listed')
end

%% Interpolated values
if M < 0.035, M = 0.035; end % Assume constant values for lower airspeeds
interp_mode = 'linear';
c_n_alpha = interp1(Mach_range,c_n_alpha_range,M,interp_mode);
S1 = interp1(Mach_range,S1_range,M,interp_mode);
S2 = interp1(Mach_range,S2_range,M,interp_mode);
alpha1_0 = interp1(Mach_range,alpha1_0_range,M,interp_mode);
alpha_ss = interp1(Mach_range,alpha_ss_range,M,interp_mode);
delta_alpha1 = interp1(Mach_range,delta_alpha1_range,M,interp_mode);
alpha_ds0 = interp1(Mach_range,alpha_ds0_range,M,interp_mode);
Df = interp1(Mach_range,Df_range,M,interp_mode);
TvL = interp1(Mach_range,TvL_range,M,interp_mode);
c_n1 = interp1(Mach_range,c_n1_range,M,interp_mode);
K0 = interp1(Mach_range,K0_range,M,interp_mode);
K1 = interp1(Mach_range,K1_range,M,interp_mode);
K2 = interp1(Mach_range,K2_range,M,interp_mode);
kappa = interp1(Mach_range,kappa_range,M,interp_mode);
Tp = interp1(Mach_range,Tp_range,M,interp_mode);
Ta = interp1(Mach_range,Ta_range,M,interp_mode);
Tf0 = interp1(Mach_range,Tf0_range,M,interp_mode);
Tv0 = interp1(Mach_range,Tv0_range,M,interp_mode);
q0 = interp1(Mach_range,q0_range,M,interp_mode);
c_m0 = interp1(Mach_range,c_m0_range,M,interp_mode);
c_d0 = interp1(Mach_range,c_d0_range,M,interp_mode);
eta = interp1(Mach_range,eta_range,M,interp_mode);
E0 = interp1(Mach_range,E0_range,M,interp_mode);
x_ac = 0.25-K0;

%% Time delay constants adjustment for dimensional time
Tp = Tp*b/U;
Tf0 = Tf0*b/U;
Tv0 = Tv0*b/U;
TvL = TvL*b/U;
Ta = Ta*b/U;

end