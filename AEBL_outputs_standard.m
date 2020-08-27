function y = AEBL_outputs_standard(t,x,xdot_g,inputs_outputs,U,b,ah,M,beta,c_n1,eta,E0,Df,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,b1,b2,b3,b4,A1,A2,A3,A4,c_m0,c_n_alpha,TvL,y0)

%% Effective AoA and pitch rate
alpha_bar = atan((U*sin(x(16))+x(13)*cos(x(16))-b*(1/2+ah)*x(14))/(U*cos(x(16))-x(13)*sin(x(16))));
alpha_bar_dot = 1/(tan(alpha_bar)^2+1)*(((U*x(14)*cos(x(16))+xdot_g(13)*cos(x(16))-x(13)*x(14)*sin(x(16))-b*(1/2+ah)*xdot_g(14))*(U*cos(x(16))-x(13)*sin(x(16)))+(U*sin(x(16))+x(13)*cos(x(16))-b*(1/2+ah)*x(14))*(U*x(14)*sin(x(16))+xdot_g(13)*sin(x(16))+x(13)*x(14)*cos(x(16))))/(U*cos(x(16))-x(13)*sin(x(16)))^2);
q = 2*alpha_bar_dot*b/U;

%% Circulatory angle of attack 
alpha_E = beta^2*U/b*(A1*b1*x(1)+A2*b2*x(2));     

%% Find time of stall onset
tv0 = inputs_outputs(1,1); so_i = inputs_outputs(2,1); 
if abs(x(9))>=c_n1 && so_i<c_n1 
    tv0 = t; 
end
tau_v = t-tv0; if tau_v < 0, tau_v = 0; end

%% Calculate coefficients
%% c_n
% Circulatory unsteady
c_nf = c_n_alpha*alpha_E*((1+sqrt(x(10)))/2)^2;

% Impulsive
c_nI = -4/(M*K_a*T_I)*x(3)-1/(M*K_q*T_I)*x(4)+4/M*alpha_bar+1/M*q;

% Total 
c_n = c_nf+x(11)+c_nI;

%% c_m
% Circulatory unsteady
if x(10)>=x(12) % alpha_bar*alpha_bar_dot>=0 -  "if x(10)>=x(12)" is to remove discontinuity at alpha_bar_dot = 0, works better than "if alpha_bar*alpha_bar_dot>=0", even though this is what Leishman suggested
    c_mf = (K0 + K1*(1-x(10)) + K2*sin(pi*x(10)^kappa))*c_nf;
else
    c_mf = (K0 + K1*(1-x(12)) + K2*sin(pi*x(12)^kappa))*c_nf;
end

% Vortex-induced
c_mv = 0;
if tau_v<=2*TvL 
    c_mv = -0.25*(1-cos(pi*tau_v/TvL))*x(11);
end

% Impulsive
c_mI = A3/(M*b3*K_aM*T_I)*x(5)+A4/(M*b4*K_aM*T_I)*x(6)+7/(12*M*K_qM*T_I)*x(8)-1/M*alpha_bar-7/(12*M)*q;

% Total 
c_m = c_mf+c_mv+c_mI+c_m0;

%% c_c
if abs(x(9))<c_n1
    c_c = eta*c_n_alpha*(alpha_E)*sin(alpha_E)*(sqrt(x(10))-E0); % E0 taken from Sheng's modification
else
    csi = Df*(abs(x(9))-c_n1);
    c_c = eta*c_n_alpha*(alpha_E)*sin(alpha_E)*(sqrt(x(10))*x(10)^csi-E0); % See Leishman and Beddoes (1986) - A Generalised Model... - page 11
end

%% Outputs
y = [c_n; c_m; c_c; alpha_bar];
if length(y)>length(y0)
    y = y(1:length(y0));
elseif length(y)<length(y0)
    y(end+1:length(y0)) = nan;
end

end