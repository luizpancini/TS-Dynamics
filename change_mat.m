clc
clear all
close all

M = 0.07730;
dt_v = [1e-5; 5e-6];
for i=1:length(dt_v)
    filename = ['typsec_Fung_M_' num2str(M,'%.5f') '_Lyapunov_dt_' num2str(dt_v(i),'%.0e') '_r_none.mat'];
    load(filename);
    t_total = 0:dt_v(i):300;
    p = 1; j = 1;
    n = round(5/dt+1);
    while p(j)<length(L1_total)
        p = [p; j*n];
        j = j+1;
    end
    p(1) = []; p(end) = [];
    L1_total(p) = [];
    close all
    figure;loglog(t_total,abs(L1_total),'k');grid;drawnow
    save(filename,'mLE','L1_total','t_total','band_tol','x0_save','tf_save','inputs_outputs_save','d_save','dt','FD_method','d0','r');
end