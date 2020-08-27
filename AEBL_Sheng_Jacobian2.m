function J = AEBL_Sheng_Jacobian2(FD_method,N,t_i,x_i,xdot_i,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0,xdotRelTol,xdot_max_it,dx)

J_x = zeros(N);
J_xdot = zeros(N);
if strcmp(FD_method,'c2')
    for j=1:N
        Dx = zeros(1,N); Dx(j) = dx;
        % df/dx
        xdot_iter = 0;
        xdot_ipDx_g = xdot_i;
        xdot_imDx_g = xdot_i;
        xdot_ipDx_eps_rel = 10*xdotRelTol;
        xdot_imDx_eps_rel = 10*xdotRelTol;
        while xdot_ipDx_eps_rel > xdotRelTol || xdot_imDx_eps_rel > xdotRelTol
            [k1_pDx,~] = AEBL_dxdt_Sheng(t_i,x_i+Dx,xdot_ipDx_g,t_i,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            [k1_mDx,~] = AEBL_dxdt_Sheng(t_i,x_i-Dx,xdot_imDx_g,t_i,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            xdot_ipDx_ng = k1_pDx;
            xdot_imDx_ng = k1_mDx;
            xdot_ipDx_eps_rel = max(abs(1-xdot_ipDx_ng./xdot_ipDx_g));
            xdot_imDx_eps_rel = max(abs(1-xdot_imDx_ng./xdot_imDx_g));
            xdot_ipDx_g = xdot_ipDx_ng;
            xdot_imDx_g = xdot_imDx_ng;
            xdot_iter = xdot_iter + 1;
            if xdot_iter == xdot_max_it
                break
            end
        end
        J_x(:,j) = (xdot_ipDx_g-xdot_imDx_g)/(2*dx);
        % df/dxdot
        xdot_iter = 0;
        xdot_ipDx_g = xdot_i;
        xdot_imDx_g = xdot_i;
        xdot_ipDx_eps_rel = 10*xdotRelTol;
        xdot_imDx_eps_rel = 10*xdotRelTol;
        while xdot_ipDx_eps_rel > xdotRelTol || xdot_imDx_eps_rel > xdotRelTol
            [k1_pDx,~] = AEBL_dxdt_Sheng(t_i,x_i,xdot_ipDx_g+Dx,t_i,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            [k1_mDx,~] = AEBL_dxdt_Sheng(t_i,x_i,xdot_imDx_g-Dx,t_i,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            xdot_ipDx_ng = k1_pDx;
            xdot_imDx_ng = k1_mDx;
            xdot_ipDx_eps_rel = max(abs(1-xdot_ipDx_ng./xdot_ipDx_g));
            xdot_imDx_eps_rel = max(abs(1-xdot_imDx_ng./xdot_imDx_g));
            xdot_ipDx_g = xdot_ipDx_ng;
            xdot_imDx_g = xdot_imDx_ng;
            xdot_iter = xdot_iter + 1;
            if xdot_iter == xdot_max_it
                break
            end
        end
        J_xdot(:,j) = (xdot_ipDx_g-xdot_imDx_g)/(2*dx);
    end  
elseif strcmp(FD_method,'c4')
    for j=1:N
        Dx = zeros(1,N); Dx(j) = dx;
        xdot_iter = 0;
        xdot_ipDx_g = xdot_i;
        xdot_imDx_g = xdot_i;
        xdot_ipDx_eps_rel = 10*xdotRelTol;
        xdot_imDx_eps_rel = 10*xdotRelTol;
        xdot_ip2Dx_g = xdot_i;
        xdot_im2Dx_g = xdot_i;
        xdot_ip2Dx_eps_rel = 10*xdotRelTol;
        xdot_im2Dx_eps_rel = 10*xdotRelTol;
        while xdot_ipDx_eps_rel > xdotRelTol || xdot_imDx_eps_rel > xdotRelTol || xdot_ip2Dx_eps_rel > xdotRelTol || xdot_im2Dx_eps_rel > xdotRelTol
            % 1 point around
            [k1_pDx,~] = AEBL_dxdt_Sheng(t_i,x_i+Dx,xdot_ipDx_g,t_i,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            [k1_mDx,~] = AEBL_dxdt_Sheng(t_i,x_i-Dx,xdot_imDx_g,t_i,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            xdot_ipDx_ng = k1_pDx;
            xdot_imDx_ng = k1_mDx;
            xdot_ipDx_eps_rel = max(abs(1-xdot_ipDx_ng./xdot_ipDx_g));
            xdot_imDx_eps_rel = max(abs(1-xdot_imDx_ng./xdot_imDx_g));
            xdot_ipDx_g = xdot_ipDx_ng;
            xdot_imDx_g = xdot_imDx_ng;
            % 2 points around
            [k1_p2Dx,~] = AEBL_dxdt_Sheng(t_i,x_i+2*Dx,xdot_ip2Dx_g,t_i,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            [k1_m2Dx,~] = AEBL_dxdt_Sheng(t_i,x_i-2*Dx,xdot_im2Dx_g,t_i,inputs_outputs,A,B,U,M,b,ah,qinf,alpha1_0,S1,S2,eta,E0,Df,TvL,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,beta,Tf0,Tv0,Ta,m,I_alpha,S_alpha,K0,K1,K2,kappa,K_a,K_aM,K_q,K_qM,T_I,c_m0,Kh,Kh3,h_F,Ka,Ka3,alpha_F,Ch,Ca,plunge_dof,alpha_ss,alpha_ds0,q0);
            xdot_ip2Dx_ng = k1_p2Dx;
            xdot_im2Dx_ng = k1_m2Dx;
            xdot_ip2Dx_eps_rel = max(abs(1-xdot_ip2Dx_ng./xdot_ip2Dx_g));
            xdot_im2Dx_eps_rel = max(abs(1-xdot_im2Dx_ng./xdot_im2Dx_g));
            xdot_ip2Dx_g = xdot_ip2Dx_ng;
            xdot_im2Dx_g = xdot_im2Dx_ng;
            xdot_iter = xdot_iter + 1;
            if xdot_iter == xdot_max_it
                break
            end
        end
        J_x(:,j) = (-xdot_ip2Dx_g+8*xdot_ipDx_g-8*xdot_imDx_g+xdot_im2Dx_g)/(12*dx);
    end
else
    error('Finite-difference method not available')
end


J = (eye(N)-J_xdot)^-1*J_x;

end