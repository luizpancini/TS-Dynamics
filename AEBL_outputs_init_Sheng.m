function y = AEBL_outputs_init_Sheng(t,x,h_dot0,alpha_0,alpha_dot0,U,b,ah,M,beta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,b1,b2,b3,b4,A1,A2,A3,A4,c_m0,c_n_alpha,eta,E0,Df,alpha_ds0,alpha_ss,q0)

%% Effective AoA and pitch rate
h_dot = h_dot0;
h_2dot = 0;
alpha = alpha_0;
alpha_dot = alpha_dot0;
alpha_2dot = 0;

alpha_bar = atan((U*sin(alpha)+h_dot*cos(alpha)-b*(1/2+ah)*alpha_dot)/(U*cos(alpha)-h_dot*sin(alpha)));
alpha_bar_dot = 1/(tan(alpha_bar)^2+1)*((U*sin(alpha)+h_dot*cos(alpha)-b*(1/2+ah)*alpha_dot)/(U*cos(alpha)-h_dot*sin(alpha))^2*(U*alpha_dot*sin(alpha)+h_2dot*sin(alpha)+h_dot*alpha_dot*cos(alpha))+(U*alpha_dot*cos(alpha)+h_2dot*cos(alpha)-h_dot*alpha_dot*sin(alpha)-b*(1/2+ah)*alpha_2dot)/(U*cos(alpha)-h_dot*sin(alpha)));
q = 2*alpha_bar_dot*b/U;

%% Stall onset criterion
if abs(q)>=q0
    alpha_cr = alpha_ds0;       % Critical AoA
else
    alpha_cr = alpha_ss + (alpha_ds0-alpha_ss)*abs(q)/q0;
end
if alpha_cr < alpha_ss || alpha_bar*alpha_bar_dot<0 ; alpha_cr = alpha_ss; end

%% Airloads coefficients
%% c_n
% Unsteady circulatory
alpha_E = beta^2*U/b*(A1*b1*x(1)+A2*b2*x(2));
c_nf = c_n_alpha*alpha_E*((1+sqrt(x(10)))/2)^2;
% Impulsive
c_nI = -4/(M*K_a*T_I)*x(3)-1/(M*K_q*T_I)*x(4)+4/M*alpha_bar+1/M*q;

% Total
c_n = c_nf+x(11)+c_nI;

%% c_m
% Impulsive
c_mI = A3/(M*b3*K_aM*T_I)*x(5) + A4/(M*b4*K_aM*T_I)*x(6) + 7/(12*M*K_qM*T_I)*x(8) - 1/M*alpha_bar - 7/(12*M)*q;

% Circulatory
if x(10)>=x(12) % alpha_bar*alpha_bar_dot>=0 -  "if x(10)>=x(12)" is to remove discontinuity at alpha_bar_dot = 0, works better than "if alpha_bar*alpha_bar_dot>=0", even though this is what Leishman suggested
    c_mf = (K0 + K1*(1-x(10)) + K2*sin(pi*x(10)^kappa))*c_nf;
else
    c_mf = (K0 + K1*(1-x(12)) + K2*sin(pi*x(12)^kappa))*c_nf;
end

% Vortex-induced
c_mv = 0;

% Total
c_m = c_mf+c_mv+c_mI+c_m0;

%% c_c
if abs(x(9))<alpha_cr
    c_c = eta*c_n_alpha*alpha_E*sin(alpha_E)*(sqrt(x(10))-E0); % E0 taken from Sheng's modification
else
    csi = Df*c_n_alpha*(abs(x(9))-alpha_cr);
    if alpha_bar*alpha_bar_dot<0
        csi = csi/(1+abs(q/(4*q0)));     % My modification: Accelerate rate of return to attached flow value
    end
    c_c = eta*c_n_alpha*alpha_E*sin(alpha_E)*(sqrt(x(10))*x(10)^csi-E0); % See Leishman and Beddoes (1986) - A Generalised Model... - page 11
end

%% Outputs
y = [c_n; c_m; c_c];
end