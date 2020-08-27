function [J,F,dF_dM,H,d2F_dxdM] = all_derivatives(N,dx,x,xdot_g,M,airfoil,a_inf,rho,b,ah,m,I_alpha,S_alpha,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F)

%% Evaluate the derivatives vector without neglecting the implicit acceleration terms
% Load airspeed dependent data
[A,B,U,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,beta,Tf0,Tv0,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,alpha_ss,alpha_ds0,q0,c_n1] = load_system_data(airfoil,M,b,a_inf,rho);
% Controllable variables
xdot_abs_tol = 1e-14;
xdot_max_it = 5e3;
relax_fac_min = 2^-4;
% Initialize
xdot_eps_new = 1;
xdot_eps = 1;
relax_fac = 1;
xdot_it = 0;
% Calculate
while xdot_eps_new > xdot_abs_tol
    xdot_ng = state_derivatives_Sheng(x,xdot_g,M,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
    xdot_eps_new = max(abs((xdot_ng-xdot_g)));
    if xdot_eps_new >= xdot_eps 
        if relax_fac > relax_fac_min
            relax_fac = relax_fac/2;
        end
    end
    xdot_eps = xdot_eps_new;
    xdot_g = relax_fac*xdot_ng+(1-relax_fac)*xdot_g; 
    xdot_it = xdot_it+1;
    if xdot_it==xdot_max_it
%         warning(['xdotAbsTol not achieved at speed M = ', num2str(M) ' (xdot_eps_abs = ', num2str(xdot_eps_new), ')'])
        break
    end
end
F = xdot_g;

%% Numerical Jacobian with centered finite differences
J = zeros(N);
relax_fac = 1;
for i=1:N
    Dx = zeros(N,1); Dx(i) = dx;
    F_pDx_eps_new = 10*xdot_abs_tol; F_mDx_eps_new = 10*xdot_abs_tol; F_pDx_eps = 10*xdot_abs_tol; F_mDx_eps = 10*xdot_abs_tol;
    F_pDx_g = F; F_mDx_g = F;
    xdot_it = 0;
    while F_pDx_eps_new > xdot_abs_tol || F_mDx_eps_new > xdot_abs_tol
        F_pDx_ng = state_derivatives_Sheng(x+Dx,F_pDx_g,M,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
        F_mDx_ng = state_derivatives_Sheng(x-Dx,F_mDx_g,M,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
        F_pDx_eps_new = max(abs((F_pDx_ng-F_pDx_g)));
        F_mDx_eps_new = max(abs((F_mDx_ng-F_mDx_g)));
        if F_pDx_eps_new >= F_pDx_eps || F_mDx_eps_new >= F_mDx_eps
            if relax_fac > relax_fac_min
                relax_fac = relax_fac/2;
            end
        end
        F_pDx_eps = F_pDx_eps_new; F_mDx_eps = F_mDx_eps_new;
        F_pDx_g = relax_fac*F_pDx_ng+(1-relax_fac)*F_pDx_g;
        F_mDx_g = relax_fac*F_mDx_ng+(1-relax_fac)*F_mDx_g;
        xdot_it = xdot_it+1;
        if xdot_it==xdot_max_it
%             warning(['xdotAbsTol not achieved at speed M=', num2str(M,'%10.4f') ' (xdot_eps_abs = ', num2str(max([F_pDx_eps; F_mDx_eps])), ')'])
            break
        end
    end
    J(:,i) = (F_pDx_g-F_mDx_g)/(2*dx);
end

%% Derivatives with respect to M
M0 = M;
% Positive disturbance
M = M0+dx;
[A,B,U,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,alpha_ss,alpha_ds0,q0,c_n1] = load_system_data(airfoil,M,b,a_inf,rho);
F_pM = state_derivatives_Sheng(x,xdot_g,M,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
% Negative disturbance
M = M0-dx;
[A,B,U,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,alpha_ss,alpha_ds0,q0,c_n1] = load_system_data(airfoil,M,b,a_inf,rho);
F_mM = state_derivatives_Sheng(x,xdot_g,M,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
dF_dM = (F_pM-F_mM)/(2*dx);

%% Hessian
H = zeros(N,N^2);
for i=1:N
    Dx = zeros(N,1); Dx(i) = dx;
    for j=1:N
        Dy = zeros(N,1); Dy(j) = dx;
        F_pDx_pDy_eps_new = 10*xdot_abs_tol; F_mDx_mDy_eps_new = 10*xdot_abs_tol; F_pDx_mDy_eps_new = 10*xdot_abs_tol; F_mDx_pDy_eps_new = 10*xdot_abs_tol;
        F_pDx_pDy_eps = 10*xdot_abs_tol; F_mDx_mDy_eps = 10*xdot_abs_tol; F_pDx_mDy_eps = 10*xdot_abs_tol; F_mDx_pDy_eps = 10*xdot_abs_tol;
        F_pDx_pDy_g = F; F_mDx_mDy_g = F; F_pDx_mDy_g = F; F_mDx_pDy_g = F;
        xdot_it = 0;
        while F_pDx_pDy_eps_new > xdot_abs_tol || F_mDx_mDy_eps_new > xdot_abs_tol || F_pDx_mDy_eps_new > xdot_abs_tol || F_mDx_pDy_eps_new > xdot_abs_tol
            F_pDx_pDy_ng = state_derivatives_Sheng(x+Dx+Dy,F_pDx_pDy_g,M,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
            F_mDx_mDy_ng = state_derivatives_Sheng(x-Dx-Dy,F_mDx_mDy_g,M,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
            F_pDx_mDy_ng = state_derivatives_Sheng(x+Dx-Dy,F_pDx_mDy_g,M,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
            F_mDx_pDy_ng = state_derivatives_Sheng(x-Dx+Dy,F_mDx_pDy_g,M,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
            F_pDx_pDy_eps_new = max(abs((F_pDx_pDy_ng-F_pDx_pDy_g)));
            F_mDx_mDy_eps_new = max(abs((F_mDx_mDy_ng-F_mDx_mDy_g)));
            F_pDx_mDy_eps_new = max(abs((F_pDx_mDy_ng-F_pDx_mDy_g)));
            F_mDx_pDy_eps_new = max(abs((F_mDx_pDy_ng-F_mDx_pDy_g)));
            if F_pDx_pDy_eps_new >= F_pDx_pDy_eps || F_mDx_mDy_eps_new >= F_mDx_mDy_eps || F_pDx_mDy_eps_new >= F_pDx_mDy_eps || F_mDx_pDy_eps_new >= F_mDx_pDy_eps
                if relax_fac > relax_fac_min
                    relax_fac = relax_fac/2;
                end
            end
            F_pDx_pDy_eps = F_pDx_pDy_eps_new; F_mDx_mDy_eps = F_mDx_mDy_eps_new; F_pDx_mDy_eps = F_pDx_mDy_eps_new; F_mDx_pDy_eps = F_mDx_pDy_eps_new;
            F_pDx_pDy_g = relax_fac*F_pDx_pDy_ng+(1-relax_fac)*F_pDx_pDy_g;
            F_mDx_mDy_g = relax_fac*F_mDx_mDy_ng+(1-relax_fac)*F_mDx_mDy_g;
            F_pDx_mDy_g = relax_fac*F_pDx_mDy_ng+(1-relax_fac)*F_pDx_mDy_g;
            F_mDx_pDy_g = relax_fac*F_mDx_pDy_ng+(1-relax_fac)*F_mDx_pDy_g;
            xdot_it = xdot_it+1;
            if xdot_it==xdot_max_it
    %             warning(['xdotAbsTol not achieved at speed M = ', num2str(M) ' (xdot_eps_abs = ', num2str(max([F_pDx_pDy_eps; F_mDx_mDy_eps; F_pDx_mDy_eps; F_mDx_pDy_eps])), ')'])
                break
            end
        end
        H(:,i*j) = (F_pDx_pDy_g-F_pDx_mDy_g-F_mDx_pDy_g+F_mDx_mDy_g)/(4*dx^2);
    end
end

%% Derivatives with respect to x and M
M0 = M;
d2F_dxdM = zeros(N);
for i=1:N
    Dx = zeros(N,1); Dx(i) = dx; DM = dx;
    F_pDx_pDM_eps_new = 10*xdot_abs_tol; F_mDx_mDM_eps_new = 10*xdot_abs_tol; F_pDx_mDM_eps_new = 10*xdot_abs_tol; F_mDx_pDM_eps_new = 10*xdot_abs_tol;
    F_pDx_pDM_eps = 10*xdot_abs_tol; F_mDx_mDM_eps = 10*xdot_abs_tol; F_pDx_mDM_eps = 10*xdot_abs_tol; F_mDx_pDM_eps = 10*xdot_abs_tol;
    F_pDx_pDM_g = F; F_mDx_mDM_g = F; F_pDx_mDM_g = F; F_mDx_pDM_g = F;
    xdot_it = 0;
    while F_pDx_pDM_eps_new > xdot_abs_tol || F_mDx_mDM_eps_new > xdot_abs_tol || F_pDx_mDM_eps_new > xdot_abs_tol || F_mDx_pDM_eps_new > xdot_abs_tol
        [A,B,U,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,alpha_ss,alpha_ds0,q0,c_n1] = load_system_data(airfoil,M0+DM,b,a_inf,rho);
        F_pDx_pDM_ng = state_derivatives_Sheng(x+Dx,F_pDx_pDM_g,M0+DM,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
        F_mDx_pDM_ng = state_derivatives_Sheng(x-Dx,F_mDx_pDM_g,M0+DM,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
        [A,B,U,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,alpha_ss,alpha_ds0,q0,c_n1] = load_system_data(airfoil,M0-DM,b,a_inf,rho);
        F_mDx_mDM_ng = state_derivatives_Sheng(x-Dx,F_mDx_mDM_g,M0-DM,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
        F_pDx_mDM_ng = state_derivatives_Sheng(x+Dx,F_pDx_mDM_g,M0-DM,U,b,ah,beta,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_q,K_aM,K_qM,T_I,K0,K1,K2,kappa,alpha1_0,S1,S2,eta,Df,E0,Ta,Tf0,Tv0,TvL,alpha_ds0,alpha_ss,q0,c_m0,m,I_alpha,S_alpha,qinf,Ca,Ch,Ka,Kh,Ka3,Kh3,alpha_F,h_F);
        F_pDx_pDM_eps_new = max(abs((F_pDx_pDM_ng-F_pDx_pDM_g)));
        F_mDx_mDM_eps_new = max(abs((F_mDx_mDM_ng-F_mDx_mDM_g)));
        F_pDx_mDM_eps_new = max(abs((F_pDx_mDM_ng-F_pDx_mDM_g)));
        F_mDx_pDM_eps_new = max(abs((F_mDx_pDM_ng-F_mDx_pDM_g)));
        if F_pDx_pDM_eps_new >= F_pDx_pDM_eps || F_mDx_mDM_eps_new >= F_mDx_mDM_eps || F_pDx_mDM_eps_new >= F_pDx_mDM_eps || F_mDx_pDM_eps_new >= F_mDx_pDM_eps
            if relax_fac > relax_fac_min
                relax_fac = relax_fac/2;
            end
        end
        F_pDx_pDM_eps = F_pDx_pDM_eps_new; F_mDx_mDM_eps = F_mDx_mDM_eps_new; F_pDx_mDM_eps = F_pDx_mDM_eps_new; F_mDx_pDM_eps = F_mDx_pDM_eps_new;
        F_pDx_pDM_g = relax_fac*F_pDx_pDM_ng+(1-relax_fac)*F_pDx_pDM_g;
        F_mDx_mDM_g = relax_fac*F_mDx_mDM_ng+(1-relax_fac)*F_mDx_mDM_g;
        F_pDx_mDM_g = relax_fac*F_pDx_mDM_ng+(1-relax_fac)*F_pDx_mDM_g;
        F_mDx_pDM_g = relax_fac*F_mDx_pDM_ng+(1-relax_fac)*F_mDx_pDM_g;
        xdot_it = xdot_it+1;
        if xdot_it==xdot_max_it
%             warning(['xdotAbsTol not achieved at speed M = ', num2str(M) ' (xdot_eps_abs = ', num2str(max([F_pDx_pDM_eps; F_mDx_mDM_eps; F_pDx_mDM_eps; F_mDx_pDM_eps])), ')'])
            break
        end
    end
    d2F_dxdM(:,i) = (F_pDx_pDU_g-F_pDx_mDM_g-F_mDx_pDM_g+F_mDx_mDM_g)/(4*dx^2);
end

end