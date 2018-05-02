function [peak,peak_err,fwhm,tau]=lorentzian_peak_fit(fft,two_mode,plotty)

%How much of peak to fit with Lorentzian
percent_peak_fit=1;

if nargin<2
    two_mode=0;
end

if nargin<3
    plotty=0;
end

fft(:,1)=fft(:,1)/10^9; %put everything in units of GHz so fit is not crazy

st_point=50; %set to cut off DC spike in fft, if necessary

[max_val,peak_ind]=max(fft(st_point:end,2));
peak_loc=fft(peak_ind,1);

%normalize the fft so the scale isn't absurd
fft(:,2)=fft(:,2)/max_val;

if percent_peak_fit~=1
    neg_go=1;
    neg_ind=peak_ind;
    while neg_go
        if fft(neg_ind,2)<=percent_peak_fit
            neg_ind_final=neg_ind;
            neg_go=0;
        else
            neg_ind=neg_ind-1;
        end
    end
    
    pos_go=1;
    pos_ind=peak_ind;
    while pos_go
        if fft(pos_ind,2)<=percent_peak_fit
            pos_ind_final=pos_ind;
            pos_go=0;
        else
            pos_ind=pos_ind+1;
        end
    end
else
    neg_ind_final=1;
    pos_ind_final=length(fft(:,1));
end

ST=[1e-4 peak_loc .01 0];
LB=[0 0 0.001 0];
UB=[1 1.2 0.1 1];

OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
TYPE=fittype('(A./((x-x0)^2+(W)^2))+C','options',OPS,'coefficients',{'A','x0','W','C'});

[f0,~]=fit(fft(neg_ind_final:pos_ind_final,1),fft(neg_ind_final:pos_ind_final,2),TYPE);
error=confint(f0,0.95); %this is the 2 sigma error

peak=f0.x0*10^9; %so answers is in Hz
peak_err=((f0.x0-error(1,2))*10^9)/2; % so answer 1 sigma in Hz
fwhm=2*f0.W*10^9; %so answer is in Hz
tau=1/(pi*fwhm); %so answer is in sec

if two_mode
    new_fft_amp=fft(:,2)-f0(fft(:,1));
    [~,peak_ind_2]=max(new_fft_amp(st_point:end));
    peak_loc_2=fft(peak_ind_2,1);
    
    ST1=[1e-4 peak_loc_2 .01 0];
    OPS1=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST1);
    TYPE1=fittype('(A./((x-x0)^2+(W)^2))+C','options',OPS1,'coefficients',{'A','x0','W','C'});
    
    [f1,~]=fit(fft(:,1),new_fft_amp,TYPE1);
    error2=confint(f1,0.95);
    
    peak2=f1.x0*10^9; %so answers is in Hz
    peak_err2=((f1.x0-error2(1,2))*10^9)/2; % so answer 1 sigma in Hz
    fwhm2=2*f1.W*10^9; %so answer is in Hz
    tau2=1/(pi*fwhm2); %so answer is in sec
    
    peak=[peak peak2];
    peak_err=[peak_err peak_err2];
    fwhm=[fwhm fwhm2];
    tau=[tau tau2];
end

if plotty
    figure()
    plot(fft(:,1),fft(:,2),'k-','LineWidth',1.25);
    hold on
    plot(fft(:,1),f0(fft(:,1)),'r--','LineWidth',1.25);
    hold on
    if two_mode
        plot(fft(:,1),f1(fft(:,1)),'b--','LineWidth',1.25);
        hold on
    end
    xlim([0 1.2])
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
    xlabel({'Frequency [GHz]'},...
        'FontUnits','points',...
        'FontSize',20,...
        'FontName','Helvetica')
end

end