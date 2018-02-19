%script to process all Ge phase data using the fixed phase method and plot the results versus the BTE predictions.
%This code used to genreate Figure 5 in the manuscript

%%% !!! Need to provide correct user directories for data file locations throughout the code or else it will not work !!! %%%
%%% There are FOUR places that this needs to be corrected

%This is the fixed null-point start position. Set to 2 for best fit based on paper results
phase_position=2;

str_mid_98='09.80um-spot';
str_mid_85='08.50um-spot_';
str_mid_64='06.40um-spot';
str_mid_55='05.50um-spot';
str_mid_48='04.80um-spot';
str_mid_36='03.60um-spot';
str_end_pos='-POS-1.txt';
str_end_neg='-NEG-1.txt';

num_meas=10; %for 3.6, 4.8, and 6.4 um
num_meas_big=30; %for 5.5, 8.5, and 9.8 um
grat98=9.7956; %Claibrated grating spacing from 2017-08-29
grat85=8.4764; %Calibrated grating spacing from 2017-08-16
grat64=6.3918; %Calibrated grating spacing from 2017-08-10
grat55=5.4891; %Calibrated grating spacing from 2017-09-04
grat48=4.7985; %Calibrated grating spacing from 2017-08-10
grat36=3.5938; %Calibrated grating spacing from 2017-08-10

grat_vec=[grat36 grat48 grat55 grat64 grat85 grat98];

%Correct for user directory
cd('.../100_001_vac_phase/');
str_base='Ge_thermal-100_001_vac_phase-';

diff_98_phase=zeros(1,num_meas_big);
diff_85_phase=zeros(1,num_meas_big);
diff_64_phase=zeros(1,num_meas);
diff_55_phase=zeros(1,num_meas_big);
diff_48_phase=zeros(1,num_meas);
diff_36_phase=zeros(1,num_meas);
diff_err_98=zeros(2,num_meas_big);
diff_err_85=zeros(2,num_meas_big);
diff_err_64=zeros(2,num_meas);
diff_err_55=zeros(2,num_meas_big);
diff_err_48=zeros(2,num_meas);
diff_err_36=zeros(2,num_meas);

for jj=1:num_meas

    if jj<10
        num_str=strcat('0',num2str(jj));
    else
        num_str=num2str(jj);
    end
    
    pos_str_64=strcat(str_base,str_mid_64,num_str,str_end_pos);
    pos_str_48=strcat(str_base,str_mid_48,num_str,str_end_pos);
    pos_str_36=strcat(str_base,str_mid_36,num_str,str_end_pos);
    neg_str_64=strcat(str_base,str_mid_64,num_str,str_end_neg);
    neg_str_48=strcat(str_base,str_mid_48,num_str,str_end_neg);
    neg_str_36=strcat(str_base,str_mid_36,num_str,str_end_neg);
    
    [diff_64_phase(jj),diff_err_64(:,jj)]=thermal_phase(pos_str_64,neg_str_64,grat64,phase_position);
    [diff_48_phase(jj),diff_err_48(:,jj)]=thermal_phase(pos_str_48,neg_str_48,grat48,phase_position);
    [diff_36_phase(jj),diff_err_36(:,jj)]=thermal_phase(pos_str_36,neg_str_36,grat36,phase_position);
end

display('Block 1 finished')

for jj=1:num_meas_big
    if jj<10
        num_str=strcat('0',num2str(jj));
    else
        num_str=num2str(jj);
    end
    
    pos_str_55=strcat(str_base,str_mid_55,num_str,str_end_pos);
    pos_str_85=strcat(str_base,str_mid_85,num_str,str_end_pos);
    pos_str_98=strcat(str_base,str_mid_98,num_str,str_end_pos);
    neg_str_55=strcat(str_base,str_mid_55,num_str,str_end_neg);
    neg_str_85=strcat(str_base,str_mid_85,num_str,str_end_neg);
    neg_str_98=strcat(str_base,str_mid_98,num_str,str_end_neg);
    
    [diff_55_phase(jj),diff_err_55(:,jj)]=thermal_phase(pos_str_55,neg_str_55,grat55,phase_position);
    [diff_85_phase(jj),diff_err_85(:,jj)]=thermal_phase(pos_str_85,neg_str_85,grat85,phase_position);
    [diff_98_phase(jj),diff_err_98(:,jj)]=thermal_phase(pos_str_98,neg_str_98,grat98,phase_position);
end

display('Block 2 (big) finished')

diff_36_001=diff_36_phase;
diff_48_001=diff_48_phase;
diff_55_001=diff_55_phase;
diff_64_001=diff_64_phase;
diff_85_001=diff_85_phase;
diff_98_001=diff_98_phase;

diff_err_36_001=diff_err_36;
diff_err_48_001=diff_err_48;
diff_err_55_001=diff_err_55;
diff_err_64_001=diff_err_64;
diff_err_85_001=diff_err_85;
diff_err_98_001=diff_err_98;

diff_ave_98_001=mean(diff_98_001);
diff_ave_85_001=mean(diff_85_001);
diff_ave_64_001=mean(diff_64_001);
diff_ave_55_001=mean(diff_55_001);
diff_ave_48_001=mean(diff_48_001);
diff_ave_36_001=mean(diff_36_001);

diff_std_98_001=std(diff_98_001);
diff_std_85_001=std(diff_85_001);
diff_std_64_001=std(diff_64_001);
diff_std_55_001=std(diff_55_001);
diff_std_48_001=std(diff_48_001);
diff_std_36_001=std(diff_36_001);

%Correct for user directory
cd('.../111_011_vac_phase/');
str_base='Ge_thermal-111_011_vac_phase-';

diff_64_phase=zeros(1,num_meas);
diff_48_phase=zeros(1,num_meas);
diff_36_phase=zeros(1,num_meas);
diff_err_64=zeros(2,num_meas);
diff_err_48=zeros(2,num_meas);
diff_err_36=zeros(2,num_meas);

for jj=1:num_meas
    if jj<10
        num_str=strcat('0',num2str(jj));
    else
        num_str=num2str(jj);
    end
    
    pos_str_64=strcat(str_base,str_mid_64,num_str,str_end_pos);
    pos_str_48=strcat(str_base,str_mid_48,num_str,str_end_pos);
    pos_str_36=strcat(str_base,str_mid_36,num_str,str_end_pos);
    neg_str_64=strcat(str_base,str_mid_64,num_str,str_end_neg);
    neg_str_48=strcat(str_base,str_mid_48,num_str,str_end_neg);
    neg_str_36=strcat(str_base,str_mid_36,num_str,str_end_neg);
    
    [diff_64_phase(jj),diff_err_64(:,jj)]=thermal_phase(pos_str_64,neg_str_64,grat64,phase_position);
    [diff_48_phase(jj),diff_err_48(:,jj)]=thermal_phase(pos_str_48,neg_str_48,grat48,phase_position);
    [diff_36_phase(jj),diff_err_36(:,jj)]=thermal_phase(pos_str_36,neg_str_36,grat36,phase_position);
end

display('Block 3 finished')

diff_36_011=diff_36_phase;
diff_48_011=diff_48_phase;
diff_64_011=diff_64_phase;

diff_err_36_011=diff_err_36;
diff_err_48_011=diff_err_48;
diff_err_64_011=diff_err_64;

diff_ave_64_011=mean(diff_64_011);
diff_ave_48_011=mean(diff_48_011);
diff_ave_36_011=mean(diff_36_011);

diff_std_64_011=std(diff_64_011);
diff_std_48_011=std(diff_48_011);
diff_std_36_011=std(diff_36_011);

%Correct for user directory
cd('.../110_111_vac_phase/');
str_base='Ge_thermal-110_111_vac_phase-';

diff_64_phase=zeros(1,num_meas);
diff_48_phase=zeros(1,num_meas);
diff_36_phase=zeros(1,num_meas);
diff_err_64=zeros(2,num_meas);
diff_err_48=zeros(2,num_meas);
diff_err_36=zeros(2,num_meas);

for jj=1:num_meas
    if jj<10
        num_str=strcat('0',num2str(jj));
    else
        num_str=num2str(jj);
    end
    
    pos_str_64=strcat(str_base,str_mid_64,num_str,str_end_pos);
    pos_str_48=strcat(str_base,str_mid_48,num_str,str_end_pos);
    pos_str_36=strcat(str_base,str_mid_36,num_str,str_end_pos);
    neg_str_64=strcat(str_base,str_mid_64,num_str,str_end_neg);
    neg_str_48=strcat(str_base,str_mid_48,num_str,str_end_neg);
    neg_str_36=strcat(str_base,str_mid_36,num_str,str_end_neg);
    
    [diff_64_phase(jj),diff_err_64(:,jj)]=thermal_phase(pos_str_64,neg_str_64,grat64,phase_position);
    [diff_48_phase(jj),diff_err_48(:,jj)]=thermal_phase(pos_str_48,neg_str_48,grat48,phase_position);
    [diff_36_phase(jj),diff_err_36(:,jj)]=thermal_phase(pos_str_36,neg_str_36,grat36,phase_position);
end

display('Block 4 finished')

diff_36_111=diff_36_phase;
diff_48_111=diff_48_phase;
diff_64_111=diff_64_phase;

diff_err_36_111=diff_err_36;
diff_err_48_111=diff_err_48;
diff_err_64_111=diff_err_64;

diff_ave_64_111=mean(diff_64_111);
diff_ave_48_111=mean(diff_48_111);
diff_ave_36_111=mean(diff_36_111);

diff_std_64_111=std(diff_64_111);
diff_std_48_111=std(diff_48_111);
diff_std_36_111=std(diff_36_111);

%%% Averaging all of them now

diff_ave_36=mean([diff_36_001 diff_36_011 diff_36_111])*10^4; %to put in cm^2/s
diff_ave_48=mean([diff_48_001 diff_48_011 diff_48_111])*10^4;
diff_ave_55=mean(diff_55_001)*10^4;
diff_ave_64=mean([diff_64_001 diff_64_011 diff_64_111])*10^4;
diff_ave_85=mean(diff_85_001)*10^4;
diff_ave_98=mean(diff_98_001)*10^4;

diff_std_36=std([diff_36_001 diff_36_011 diff_36_111])*10^4;
diff_std_48=std([diff_48_001 diff_48_011 diff_48_111])*10^4;
diff_std_55=std(diff_55_001)*10^4;
diff_std_64=std([diff_64_001 diff_64_011 diff_64_111])*10^4;
diff_std_85=std(diff_85_001)*10^4;
diff_std_98=std(diff_98_001)*10^4;

%Correct for user directory
cd('.../BTE_variational_k/');

diff_vec=[diff_ave_36 diff_ave_48 diff_ave_55 diff_ave_64 diff_ave_85 diff_ave_98];
diff_err_vec=[diff_std_36 diff_std_48 diff_std_55 diff_std_64 diff_std_85 diff_std_98];

get_Ge_k(0,1);
hold on
errorbar(grat_vec,diff_vec,diff_err_vec,'k.','LineWidth',1.25)
hold on
plot([grat_vec(3) grat_vec(5:6)],[diff_ave_55 diff_ave_85 diff_ave_98],'ro','MarkerSize',8,'MarkerFaceColor','r')
hold on
plot([grat_vec(1:2) grat_vec(4)],[diff_ave_36 diff_ave_48 diff_ave_64],'rd','MarkerSize',8,'MarkerFaceColor','r')
hold on
set(gca,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',16,...
    'FontName','Helvetica',...
    'LineWidth',1.25)
ylabel('Effective Thermal Diffusivity [cm^2/s]',...
    'FontUnits','points',...
    'FontName','Helvetica',...
    'FontSize',20)
xlabel(['Grating spacing [' 956 'm]'],...
    'FontUnits','points',...
    'FontName','Helvetica',...
    'FontSize',20)
