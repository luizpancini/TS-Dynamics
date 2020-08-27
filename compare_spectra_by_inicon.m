clc
clear all
close all

typsec_data_inicon1 = 'typsec_Fung_inicon1.mat';
typsec_data_inicon2 = 'typsec_Fung_inicon2.mat';

load(typsec_data_inicon1,'f_vec','b','U','SP1_xi','SP1_alpha');
SP1_xi1 = SP1_xi; SP1_alpha1 = SP1_alpha;
load(typsec_data_inicon2,'f_vec','b','U','SP1_xi','SP1_alpha');
SP1_xi2 = SP1_xi; SP1_alpha2 = SP1_alpha;

SP1_xi_dif = abs(SP1_xi2-SP1_xi1);
SP1_alpha_dif = abs(SP1_alpha2-SP1_alpha1);

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'FontSize',16,'FontName','times new roman','YScale','log');
hold(axes1,'on');
semilogy(f_vec*2*pi*b/U,SP1_xi_dif,'k-','Parent',axes1);
xlabel('$k$','LineWidth',1,'FontWeight','normal','FontSize',16);
ylabel('Difference in amplitude spectrum of $\xi$','LineWidth',1,'FontWeight','normal','FontSize',16);
grid off

figure2 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes2 = axes('Parent',figure2,'FontSize',16,'FontName','times new roman','YScale','log');
hold(axes2,'on');
semilogy(f_vec*2*pi*b/U,SP1_alpha_dif,'k-','Parent',axes2);
xlabel('$k$','LineWidth',1,'FontWeight','normal','FontSize',16);
ylabel('Difference in amplitude spectrum of $\alpha$','LineWidth',1,'FontWeight','normal','FontSize',16);
grid off