clc
clear all
close all

M_vec = 0.07651:1e-5:0.07654;
for i=1:length(M_vec)
    M = M_vec(i);
    filename = ['typsec_Fung_M_' num2str(M,'%.5f') '_spectrum_and_poincare.mat'];
    load(filename);
    save (filename,'M','dt','f_vec','SP1_alpha','SP1_xi','SP1_alpha_bar','alpha_P','alpha_dot_P','xi_P','xi_dot_P','alpha_bar_P','alpha_bar_dot_P'); 
end