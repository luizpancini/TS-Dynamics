function [boundary,inputs_outputs] = BL_boundaries_Floquet(t_ip1,t_i,x9_i,x10_ip1,x10_i,x12_ip1,x12_i,inputs_outputs,mods,alpha1_0,c_n_alpha,TvL,q0,c_n1,xdot_i,xdot_ip1,marker)

tv0 = inputs_outputs(1,1);
inputs_outputs(2,1) = abs(x9_i); so_i = abs(x9_i);
so_ip1 = inputs_outputs(2,2); 
alpha1_i = inputs_outputs(3,1);
alpha1_ip1 = inputs_outputs(3,2);
alpha_bar_i = inputs_outputs(4,1);
alpha_bar_ip1 = inputs_outputs(4,2);
alpha_bar_dot_i = inputs_outputs(5,1);
alpha_bar_dot_ip1 = inputs_outputs(5,2);
q_i = inputs_outputs(6,1);
q_ip1 = inputs_outputs(6,2);
so_lim_i = inputs_outputs(7,1);
so_lim_ip1 = inputs_outputs(7,2);

% Discontinuity boundaries
boundary(1)=(abs(alpha_bar_ip1)-alpha1_0)*(abs(alpha_bar_i)-alpha1_0);        % f: alpha_bar and alpha1_0
if strcmp(mods,'off') % Standard model
    boundary(2)=(so_ip1-c_n1)*(so_i-c_n1);                                    % so_state and so_lim 
    boundary(3)=(so_ip1/c_n_alpha-alpha1_ip1)*(so_i/c_n_alpha-alpha1_i);      % f': x9/c_n_alpha and alpha_1
    breakpoint_sep = 0.7;
else
    boundary(2)=(so_ip1-so_lim_ip1)*(so_i-so_lim_i);                              % so_state and so_lim
    boundary(3)=(so_ip1-alpha1_ip1)*(so_i-alpha1_i);                              % f': x9 and alpha_1
    breakpoint_sep = 0.7;
    boundary(9) = (abs(q_ip1)-q0)*(abs(q_i)-q0);
    boundary(10) = q_ip1*q_i;
end
boundary(4)=(t_ip1-tv0-TvL)*(t_i-tv0-TvL);
boundary(5)=(t_ip1-tv0-2*TvL)*(t_i-tv0-2*TvL);
boundary(6)=(alpha_bar_ip1*alpha_bar_dot_ip1)*(alpha_bar_i*alpha_bar_dot_i);
boundary(7)=(x10_ip1-x12_ip1)*(x10_i-x12_i);
boundary(8)=(x10_ip1-breakpoint_sep)*(x10_i-breakpoint_sep);
if marker > 0
    boundary(11) = xdot_i(marker)*xdot_ip1(marker);
else
    boundary(11) = 1;
end
end