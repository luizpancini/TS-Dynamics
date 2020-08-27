function [df_dx,df_dT] = df_dx_and_df_dT(N,T0,x0,xdot0,f,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options,diffs)

dx = 1e-10;
dt = dx;
% Apply finite-differences to calculate df_dx
df_dx = zeros(N);
for i=1:N
    Dx = zeros(1,N); Dx(i) = dx;
    [~,x_pDx] = AEBL_rkf_Sheng_Shooting([0 T0],x0+Dx,xdot0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
    f_pDx = (x_pDx(end,:)-x_pDx(1,:))';
    if strcmp(diffs,'f')
        df_dx(:,i) = (f_pDx-f)/dx;
    elseif strcmp(diffs,'c')
        [~,x_mDx] = AEBL_rkf_Sheng_Shooting([0 T0],x0-Dx,xdot0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
        f_mDx = (x_mDx(end,:)-x_mDx(1,:))';
        df_dx(:,i) = (f_pDx-f_mDx)/(2*dx);
    end
end

% Apply centered-differences to calculate df_dT
[~,x_pDt] = AEBL_rkf_Sheng_Shooting([0 T0+dt],x0,xdot0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
[~,x_mDt] = AEBL_rkf_Sheng_Shooting([0 T0-dt],x0,xdot0,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,mods,c_n1,options);
f_pDt = (x_pDt(end,:)-x_pDt(1,:))';
f_mDt = (x_mDt(end,:)-x_mDt(1,:))';
df_dT = (f_pDt-f_mDt)/(2*dt);