function []=thermal_model_plot(grat,alpha,beta)
%Take a value of alpha and beta in, and plot the pure thermal model for
%both displacement and refectivity as well as the composite from 0 to 50 ns.

t=(0.1:.25:50)*10^(-9);
q=2*pi/(grat*10^(-6));

displacement=erfc(q*sqrt(alpha.*t));
reflectivity=beta*exp(-q^2*alpha.*t)./sqrt(t);
total=displacement-reflectivity;
[~,t0_ind]=max(total);
t0=t(t0_ind)/1e-9;

figure()
plot([t0 t0],[-0.25 1],'k--','LineWidth',1.25)
hold on
plot(t/1e-9,displacement,'b-','LineWidth',1.6)
hold on
plot(t/1e-9,reflectivity,'r-','LineWidth',1.6)
hold on
plot(t/1e-9,total,'k-','LineWidth',1.6)
hold on
ylim([-0.25 1])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Helvetica',...
    'LineWidth',1.25)
ylabel({'Amplutide [a.u.]'},...
    'FontUnits','points',...
    'FontName','Helvetica',...
    'FontSize',20)
xlabel({'Time [ns]'},...
    'FontUnits','points',...
    'FontName','Helvetica',...
    'FontSize',20)

end