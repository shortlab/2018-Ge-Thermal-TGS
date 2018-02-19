%read and plot the start variation 

data_without=dlmread('start_var_comp_2nd_6.txt');
data_with=dlmread('sine_start_var_comp_2nd_6.txt');

start_time=2.05; %in ns
period=2.25; %in ns

start_index=start_time/0.05-1;
period_index=period/0.05;
ave_start_index=floor(start_index-floor(period_index/2));
ave_end_index=floor(start_index+ceil(period_index/2));

data_without(:,1)=data_without(:,1)*10^9; %to put in units of ns
data_with(:,1)=data_with(:,1)*10^9; %to put in units of ns
data_without(:,2)=data_without(:,2)*10^4; %to put in units of cm^2/s
data_with(:,2)=data_with(:,2)*10^4; %to put in units of cm^2/s

ave_without=mean(data_without(ave_start_index:ave_end_index,2));
std_without=std(data_without(ave_start_index:ave_end_index,2));
ave_with=mean(data_with(ave_start_index:ave_end_index,2));
std_with=std(data_with(ave_start_index:ave_end_index,2));

figure()
plot([start_time start_time],[0.18 0.32],'r-.','LineWidth',1.1);
hold on
plot(data_without(:,1),data_without(:,2),'k-','LineWidth',1.25);
hold on
plot(data_with(:,1),data_with(:,2),'k--','LineWidth',1.25);
hold on
ylim([0.18 0.32])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Helvetica',...
    'LineWidth',1.25)
ylabel({'Effective Thermal  Diffusivity [cm^2/s]'},...
    'FontUnits','points',...
    'FontName','Helvetica',...
    'FontSize',20)
xlabel({'Time [ns]'},...
    'FontUnits','points',...
    'FontName','Helvetica',...
    'FontSize',20)

display(ave_without)
display(std_without/ave_without)
display(ave_with)
display(std_with/ave_with)