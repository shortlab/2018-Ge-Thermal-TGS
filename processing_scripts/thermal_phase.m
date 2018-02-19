function [diffusivity,diffusivity_err] = thermal_phase(pos_file,neg_file,grat,start_phase,end_time)
%   Function to determine thermal diffusivity from phase grating TGS data
%   Data is saved in two files, positive (with one heterodyne phsae) and negative (with another), must provide both files
%   pos_file: positive phase TGS data file
%   neg_file: negative phae TGS data file
%   grat: calibrated grating spacing in um
%   start_phase: provide integer between 1 and 4 to pick the null-point start from which fit will begin
%   end_time: a shortened fit end time if you do not want ot fit the whole profile. 
%               If argument not given, will be set to default for 200ns data collection window

%%%%Write this to include a sine variation in the fit by default, but to
%%%%start the fits from a fixed null point, not time, relative to
%%%%the initial SAW maximum. 

%Settings for various plotting and output options to be set by boolean arguments
find_max=0;
plotty=0;
plot_trace=0;
plot_final=0;
print_final_fit=0;
two_detectors=1;
q=2*pi/(grat*10^(-6));
tstep=5e-11; %Set by scope used for data collection
no_pre_calc=0;

derivative=0;

%How far from guessed values for diffusivity and beta do you vary in the
%end, d for diffusivity b for beta. Diffusivity is most important
percent_range_d=0.45;
percent_range_b=2.1;

if nargin<5
    %     end_time=10e-7; %for 50ns base on scope
    end_time=2e-7; %for 20ns base on scope
end

%Difference in file write format based on newer or older acquisition. hdr_len should be 16 for the Ge dataset
if two_detectors
    hdr_len=16;
else
    hdr_len=15;
end

%Generate filtered power spectrum using make_fft_embed_time and find the peak frequency from that profile using fit_spectra_peaks
[fft,~]=make_fft_embed_time(1,pos_file,neg_file,'therm',5,grat,0);
[peak_freq,~,~]=fit_spectra_peaks(fft,0);

%read in data files for this procedure
pos=dlmread(pos_file,'',hdr_len,0);
neg=dlmread(neg_file,'',hdr_len,0);

%sometimes written data is off by one time step at the end, chop that off if they do not match
if length(pos(:,1))>length(neg(:,1))
    pos=pos(1:length(neg(:,1)),:);
elseif length(neg(:,1))>length(neg(:,1))
    neg=neg(1:length(pos(:,1)),:);
end

%normalize each set of data to the zero level before the pump impulse
pos(:,2)=pos(:,2)-mean(pos(1:50,2));
neg(:,2)=neg(:,2)-mean(neg(1:50,2));

time_index=186; %From peak in amp grating data
time_naught=neg(time_index,1);

end_index=floor(end_time/tstep);

%re-normalize data to end signal decayed state if grating decays entirely during collection window, if not, do not re-normalize
if grat<8
    base_index=floor(length(pos(:,2))/5);
    long_base=mean((pos(end-base_index:end,2)-neg(end-base_index:end,2)));
else
    long_base=0;
end

if plot_trace
    figure()
    plot(neg(:,1)*10^9,(pos(:,2)-neg(:,2)-long_base)*10^3,'-','Color',[0 0 0.75],'LineWidth',1.25)
    hold on
    plot([neg(1,1) neg(end,1)]*10^9,[0 0],'k--','LineWidth',1.5)
    xlim([0 end_time*10^9])
    set(gca,...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',16,...
        'FontName','Helvetica',...
        'LineWidth',1.25)
    ylabel({'Amplitude [mV]'},...
        'FontUnits','points',...
        'FontSize',20,...
        'FontName','Helvetica')
    xlabel({'Time [ns]'},...
        'FontUnits','points',...
        'FontSize',20,...
        'FontName','Helvetica')
end

fixed_short=[pos(time_index:end_index,1)-time_naught pos(time_index:end_index,2)-neg(time_index:end_index,2)-long_base];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%Plot block for fixed null point figure (Figure 4)
% null_pts=4;
% null_pt_vec_x=zeros(1,null_pts);
% null_pt_vec_y=zeros(1,null_pts);
% for kk=1:length(null_pt_vec_x)
%     temp=find_start_phase(fixed_short(:,1),fixed_short(:,2),kk,grat);
%     temp_ind=floor(temp/tstep)+1;
%     null_pt_vec_x(kk)=fixed_short(temp_ind,1);
%     null_pt_vec_y(kk)=fixed_short(temp_ind,2);
% end
% 
% null_pt_vec_x=null_pt_vec_x*10^9; %scale for units
% null_pt_vec_y=null_pt_vec_y*10^3;
% 
% display(null_pt_vec_x)
% 
% figure()
% plot(fixed_short(:,1)*10^9,fixed_short(:,2)*10^3,'-','Color',[0 0 0.75],'LineWidth',1.25)
% hold on
% plot(null_pt_vec_x,null_pt_vec_y,'rd','MarkerSize',8,'MarkerFaceColor','r')
% hold on
% xlim([0 5])
% ylim([0 65])
% set(gca,...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'FontSize',16,...
%         'FontName','Helvetica',...
%         'LineWidth',1.25)
%     ylabel({'Amplitude [mV]'},...
%         'FontUnits','points',...
%         'FontSize',20,...
%         'FontName','Helvetica')
%     xlabel({'Time [ns]'},...
%         'FontUnits','points',...
%         'FontSize',20,...
%         'FontName','Helvetica')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if derivative
    der_len=length(fixed_short(:,1))-1;
    fixed_derivative=zeros(1,der_len);
    for jj=1:der_len
        fixed_derivative(jj)=(fixed_short(jj+1,2)-fixed_short(jj,2))/tstep;
    end
    figure()
    plot(fixed_short(1:der_len,1),fixed_derivative,'k-')
    title('This is the derivative of fixed short')
end

%if you don't want to automatically find the peak t_0 of the profile, setting find_max to true above will allow
%you to select a region on the plot within with to search. Useful if there are initial transients.
if find_max || plotty
    figure()
    plot(fixed_short(:,1),fixed_short(:,2),'k-')
    xlim([0 1.5*10^-7])
    title('this is fixed short');
    if find_max
        hold on
        [x_cord,~]=ginput(2);
        neg_x_cord=x_cord(1);
        pos_x_cord=x_cord(2);
        pos_x_ind=floor(pos_x_cord/tstep);
        neg_x_ind=floor(neg_x_cord/tstep);
        [~,max_index]=max(fixed_short(neg_x_ind:pos_x_ind,2));
        time_max=fixed_short(max_index+neg_x_ind,1);
        close(gcf)
    end
end

%Otherwise, find t_0 from the profile directly
if ~find_max
    time_offset_index=20; %was 20 before
    [~,max_index]=max(fixed_short(time_offset_index:end,2));
    max_index=max_index+time_offset_index-1;
    time_max=fixed_short(max_index,1);
end

start_time_master=find_start_phase(fixed_short(:,1),fixed_short(:,2),start_phase,grat);
start_index_master=floor(start_time_master/tstep)+1;

%%%%%%%%%
% %Block to write the start time for each measurement based on the fixed null
% cd('start_times')
% dlmwrite(strcat('start_time_null_2_',num2str(grat),'.txt'),start_time_master,'-append')
% cd('..')
%%%%%%%%%

%Fitting parameters for initial naive fit
LB=[0 0];
UB=[1 10^-4];
ST=[.05 5*10^-6];

OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
TYPE=fittype('A.*erfc(q*sqrt(k*(x+time_max)))','options',OPS,'problem',{'q','time_max'},'coefficients',{'A','k'});
[f0,gof]=fit(fixed_short(:,1),fixed_short(:,2),TYPE,'problem',{q,time_max});

diffusivity=f0.k;
error=confint(f0,0.95);
diffusivity_err=[diffusivity-error(1,2) error(2,2)-diffusivity];

if plotty
    figure()
    plot(fixed_short(:,1),fixed_short(:,2),fixed_short(:,1),f0(fixed_short(:,1)))
    hold on
    title('First naive fit')
end

%We'll call the parameter beta the ratio of the amplitudes of the
%displacement versus temperature grating. beta should be a small number.

for jj=1:10
    
    beta=q*sqrt(diffusivity/pi)*(q^2*diffusivity+1/(2*time_max))^(-1);
    
    start_time=start_time_master;
    start_index=start_index_master;

    %Conduct initial parameter estimation without using an sin(x) contribution to the fit
    
    LB1=[0 0];
    UB1=[1 10^-4];
    ST1=[.05 10^-5];
    
    OPS1=fitoptions('Method','NonLinearLeastSquares','Lower',LB1,'Upper',UB1,'Start',ST1);
    TYPE1=fittype('A.*(erfc(q*sqrt(k*(x+start_time)))-beta*exp(-q^2*k*(x+start_time))./sqrt((x+start_time)))','options',OPS1,'problem',{'q','beta','start_time'},'coefficients',{'A','k'});
    [f1,gof]=fit(fixed_short(start_index:end,1),fixed_short(start_index:end,2),TYPE1,'problem',{q,beta,start_time});
    
    diffusivity=f1.k;
    error=confint(f1,0.95);
    diffusivity_err=[diffusivity-error(1,2) error(2,2)-diffusivity];
    
    if plotty
        figure()
        plot(fixed_short(:,1),fixed_short(:,2),fixed_short_old(:,1),f1(fixed_short_old(:,1)))
        hold on
        title(strcat('Fit number ',num2str(jj+1),' - fixed beta'))
    end
    
end

%If you've elected not to pre-compute, provide hard initial guesses for diffusivity and beta based on 
%the bulk value of Ge diffusivity. Set ranges for final fit.
if no_pre_calc
    diffusivity=0.3636*10^-4;
    beta=2e-5;
    low_bound=[1e-5 0];
    up_bound=[1e-3 1e-4];
else
    low_bound=[diffusivity*(1-percent_range_d) beta*(1-percent_range_b)];
    up_bound=[diffusivity*(1+percent_range_d) beta*(1+percent_range_b)];
    if percent_range_d>1
        low_bound(1)=0;
    end
    if percent_range_b>1
        low_bound(2)=0;
    end
end


start_time2=start_time_master;
start_index2=start_index_master;

%After testing, initial guesses for acoustic damping are more appropriate for different grating spacings
LB2=[0 low_bound(1) low_bound(2) 0 -2*pi 0 -5e-3];
if grat<4
    UB2=[1 up_bound(1) up_bound(2) 10 2*pi 2e-7 5e-3];
    ST2=[.05 diffusivity beta 0.05 0 15e-8 0];
else
    UB2=[1 up_bound(1) up_bound(2) 10 2*pi 1e-7 5e-3];
    ST2=[.05 diffusivity beta 0.05 0 1e-8 0];
end

OPS2=fitoptions('Method','NonLinearLeastSquares','Lower',LB2,'Upper',UB2,'Start',ST2);
TYPE2=fittype('A.*(erfc(q*sqrt(k*(x+start_time)))-beta*exp(-q^2*k*(x+start_time))./sqrt((x+start_time)))+B.*sin(2*pi*(peak_freq)*(x+start_time)+p)*exp(-(x+start_time)/t)+D','options',OPS2,'problem',{'q','start_time','peak_freq'},'coefficients',{'A','k','beta','B','p','t','D'});
[f2,gof]=fit(fixed_short(start_index2:end,1),fixed_short(start_index2:end,2),TYPE2,'problem',{q,start_time2,peak_freq});

if print_final_fit
    display(f2)
end

fit_err=fixed_short(:,2)-f2(fixed_short(:,1));
sum_squared_err=sum(fit_err.^2);
deg_free=length(fixed_short(:,1))-3;
rms_err=sqrt(sum_squared_err/deg_free);

diffusivity=f2.k;

error=confint(f2,0.95);
diffusivity_err=[diffusivity-error(1,2) error(2,2)-diffusivity];

%first checks that diffusivity has not pegged to fit bounds, second checks that beta
%has not pegged. If either do, it is a bad fit
if isnan(diffusivity_err(1))
    display(strcat('Bad fit for: ',pos_file,'~re: alpha'))
elseif isnan(error(1,3))
    display(strcat('Bad fit for: ',pos_file,'~re: beta'))
end

if plot_final
    %Plotting factor for generation of traces in Figure 6
    amp_factor=1;
    %%%%%%%%%%%%%%%
    %Block to reconstruct the best-fit model without the sinusoidal
    %contribution, for comparison
    f_remove_sine=cfit(TYPE2,f2.A,f2.k,f2.beta,f2.B,f2.p,0,f2.D,q,start_time2,peak_freq);
    %%%%%%%%%%%%%%%
    figure()
%     plot((fixed_short(:,1))*10^9,(fixed_short(:,2))*10^3,'k-','LineWidth',1.5)
    plot((neg(:,1)-time_naught)*10^9,(pos(:,2)-neg(:,2)-long_base)*10^3/amp_factor,'k-','LineWidth',1.35)
    hold on
    %plot vertical line at start time
%     plot([fixed_short(start_index_master,1) fixed_short(start_index_master,1)]*10^9,ylim,'b--')
%     hold on
    plot(fixed_short(:,1)*10^9,(f2(fixed_short(:,1)))*10^3/amp_factor,'r--','LineWidth',1.45)
    hold on
    plot(fixed_short(:,1)*10^9,(f_remove_sine(fixed_short(:,1)))*10^3/amp_factor,'-','Color',[0 0 0.75],'LineWidth',1.45)
    hold on
    xlim([-5 end_time*10^9])
    set(gcf,'Position',[0 0 400 300])
    hold on
    set(gca,...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',16,...
        'FontName','Helvetica',...
        'LineWidth',1.25)
%         'YTickLabel','')
    ylabel({'Amplitude [a.u.]'},...
        'FontUnits','points',...
        'FontSize',20,...
        'FontName','Helvetica')
    xlabel({'Time [ns]'},...
        'FontUnits','points',...
        'FontSize',20,...
        'FontName','Helvetica')
end

end