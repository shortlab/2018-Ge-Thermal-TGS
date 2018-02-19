%script to plot representative phase grating and amplitude grating
%measurements

%%% !!! Must change user directory to properly locate data files in TWO places for script to run properly !!! %%%

phase_str_pos='Ge_thermal-111_011_vac_phase-04.80um-spot07-POS-1.txt';
phase_str_neg='Ge_thermal-111_011_vac_phase-04.80um-spot07-NEG-1.txt';
amp_str_pos='Ge_thermal-111_011_vac_amp-04.80um-spot01-POS-1.txt';
amp_str_neg='Ge_thermal-111_011_vac_amp-04.80um-spot01-NEG-1.txt';

hdr_len=16;
end_time=100;

%Change user directory to run correctly
cd('.../111_011_vac_phase/')

[fft,~]=make_fft_embed_time(1,phase_str_pos,phase_str_neg,'therm',5,5.5,0);

pos=dlmread(phase_str_pos,'',hdr_len,0);
neg=dlmread(phase_str_neg,'',hdr_len,0);

pos(:,2)=pos(:,2)-mean(pos(1:50,2));
neg(:,2)=neg(:,2)-mean(neg(1:50,2));

base_index=floor(length(pos(:,2))/5);
long_base=mean((pos(end-base_index:end,2)-neg(end-base_index:end,2)));

phase_x=neg(:,1)*10^9;
phase_y=(pos(:,2)-neg(:,2)-long_base)*10^3;
phase_y=phase_y*3;

%Change user directory to run properly
cd('.../111_011_vac_amp/')

[fft_a,~]=make_fft_embed_time(1,amp_str_pos,amp_str_neg,'therm',5,5.5,0);

pos_a=dlmread(amp_str_pos,'',hdr_len,0);
neg_a=dlmread(amp_str_neg,'',hdr_len,0);

pos_a(:,2)=pos_a(:,2)-mean(pos_a(1:50,2));
neg_a(:,2)=neg_a(:,2)-mean(neg_a(1:50,2));

base_index_a=floor(length(pos_a(:,2))/5);
long_base_a=mean((pos_a(end-base_index:end,2)-neg_a(end-base_index:end,2)));

amp_x=neg_a(:,1)*10^9;
amp_y=(pos_a(:,2)-neg_a(:,2)-long_base_a)*10^3;

figure()
plot([phase_x(1) phase_x(end)],[0 0],'k--','LineWidth',1.25)
hold on
plot(phase_x,phase_y,'-','Color',[0 0 0.75],'LineWidth',1.4)
hold on
plot(amp_x,amp_y,'r-','LineWidth',1.4);
hold on
xlim([0 end_time])
%     ylim([-3.5 62])
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Helvetica',...
    'LineWidth',1.25)
ylabel({'Amplitude [a.u.]'},...
    'FontUnits','points',...
    'FontSize',20,...
    'FontName','Helvetica')
xlabel({'Time [ns]'},...
    'FontUnits','points',...
    'FontSize',20,...
    'FontName','Helvetica')

figure()
fft(:,1)=fft(:,1)/10^9;
fft(:,2)=fft(:,2)/(5*10^(-8));
fft_a(:,1)=fft_a(:,1)/10^9;
fft_a(:,2)=fft_a(:,2)/(5*10^(-8));
plot(fft(:,1),fft(:,2),'-','Color',[0 0 0.75],'LineWidth',1.25)
hold on
plot(fft_a(:,1),fft_a(:,2),'r-','LineWidth',1.25)
hold on
xlim([0 1.2])
ylim([0 1]);
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Helvetica',...
    'LineWidth',1.25)
ylabel({'Power [a.u.]'},...
    'FontUnits','points',...
    'FontSize',20,...
    'FontName','Helvetica')
xlabel({'Frequency [GHz]'},...
    'FontUnits','points',...
    'FontSize',20,...
    'FontName','Helvetica')