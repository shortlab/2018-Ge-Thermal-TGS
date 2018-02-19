function [fft,diffusivity] = make_fft_embed_time(num_in,pos_file,neg_file,mod,strt,grat,peak)
%   For this use only, generate filtered fft for input data
%   num_in: legacy command, set to unity for all processing
%   pos_file: positive phase data pos_file
%   neg_file: negative phase data neg_file
%   mod: provide model with which to subtract background for filting trace before fft, read code to see options
%        only relevant choice for the purpose here is the 'therm' option, will chose this by default if none slected
%   strt: no longer used. set to 5 or 6 for processing, default to 6 if not provided
%   grat: grating, either calibrated or estimated, in um
%   peak: 

%Boolean options for plotting and processing
derivative=0;
plotty=0;
plotfft=0;
saveout=0;

if nargin==3
    mod='therm';
    strt=5;
end

if nargin==4
    strt=6;
end

hdr_len=16;

pos=dlmread(pos_file,'',hdr_len,0);
neg=dlmread(neg_file,'',hdr_len,0);

pos(:,2)=pos(:,2)-mean(pos(1:50,2));
neg(:,2)=neg(:,2)-mean(neg(1:50,2));

[~,time_index]=max(neg(1:1000,2));
time_naught=neg(time_index,1);

%If recorded traces differ in length, fix them to shorter of the two
if length(pos(:,1))>length(neg(:,1))
    pos=pos(1:length(neg(:,1)),:);
elseif length(pos(:,1))<length(neg(:,1))
    neg=neg(1:length(pos(:,1)),:);
end

fixed_short=[pos(:,1) pos(:,2)-neg(:,2)];

[~,fix_index]=max(fixed_short(:,2));
fixed_short=fixed_short(fix_index:end,:);

if saveout
    dlmwrite('dat_temp.txt',fixed_short);
end

if plotty
    figure()
    plot(fixed_short(:,1),fixed_short(:,2),'r')
    title('this is fixed_short');
end

%normalize to initial level before pump pulse
fixed_short(:,2)=fixed_short(:,2)-mean(fixed_short(end-50:end,2));

%Model option controls for trace filtering
custom=0;
therm=0;
thermsaw=0;
thermsaw_mod=0;
thermsaw2=0;

if strcmp(mod,'exp1') %Single Exponential
    mod_str=mod;
    fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
elseif strcmp(mod,'exp2') %Sum of Exponentials
    mod_str=mod;
    fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
elseif strcmp(mod,'power2') % Two term power fit
    mod_str=mod;
elseif strcmp(mod,'qexp') %Forced two term decaying exponential fit
    custom=1;
    fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
elseif strcmp(mod, 'therm') %fit thermal decay to find diffusivity
    therm=1;
    fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
elseif strcmp(mod, 'thermsaw')
    thermsaw=1;
%     thermsaw_mod=1;
    fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
elseif strcmp(mod, 'thermsaw2')
    if length(peak)<=1
        thermsaw=1;
    else
        thermsaw2=1;
    end
    fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily
else
    mod_str='power1'; %single term power fit (goes something like 1/x)
end

%Fit data to one of a couple of models determined by the mod argument in
%the function

if custom
    fixed_short(:,1)=fixed_short(:,1)-time_naught;
    q=2*pi/(grat*10^(-6));
    time_offset=fixed_short(strt,1);

    LB=[0 0];
    UB=[1 10^-4];
    ST=[.05 10^-5];
    
    OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
    TYPE=fittype('A.*erfc(q*sqrt(k*(x+time0)))','options',OPS,'problem',{'q','time0'},'coefficients',{'A','k'});

    [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',{q,time_offset});
    diffusivity=f0.k;
elseif therm
    q=2*pi/(grat*10^(-6));
    LB=[0 0 0];
    UB=[10 5*10^-4 0.1];
    ST=[0.05 5*10^-5 0];
    OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
    TYPE=fittype('A.*erfc(q*sqrt(k*(x+2.5e-9)))+c;','options',OPS,'problem','q','coefficients',{'A','k','c'});
    
    [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',q);
    diffusivity=f0.k;
elseif thermsaw
    q=2*pi/(grat*10^(-6));
    LB=[0 0 0 0 0 0];
    UB=[10 10^-3 10 2*pi 10^-6 1];
    ST=[0.005 5*10^-5 0.001 pi 10^-5 0];
    OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
    TYPE=fittype('A.*erfc(q*sqrt(k*(x+2.5e-9)))+B.*sin(2*pi*f*(x+2.5e-9)+p)*exp(-(x+2.5e-9)/t)+C;','options',OPS,'problem',{'q','f'},'coefficients',{'A','k','B','p','t','C'});    
    
    [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',{q,peak});

    diffusivity=f0.k;

elseif thermsaw_mod
    q=2*pi/(grat*10^(-6));
    LB=[0 0 0 0 0 0 0];
    UB=[10 10^-3 1 10 2*pi 10^-6 1];
    ST=[0.005 5*10^-5 0.05 0.01 0 10^-5 0];
    OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
    TYPE=fittype('A.*(erfc(q*sqrt(k*(x+2.5e-9)))-D.*exp(q^2*k*(x+2.5e-9)))+B.*sin(2*pi*f*x+p)*exp(-x/t)+C;','options',OPS,'problem',{'q','f'},'coefficients',{'A','k','D','B','p','t','C'});
    
    [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',{q,peak});
    diffusivity=f0.k;

elseif thermsaw2
    q=2*pi/(grat*10^(-6));
    LB=[0 0 0 0 0 0 0 0 0];
    UB=[10 10^-3 10 2*pi 10^-6 10 2*pi 10^-6 1];
    ST=[0.05 0.5*10^-4 0.01 0 10^-5 0.01 0 10^-5 0];
    OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
    TYPE=fittype('A.*erfc(q*sqrt(k*(x+2.5e-9)))+B.*sin(2*pi*f1*(x+2.5e-9)+p1)*exp(-(x+2.5e-9)/t1)+C.*sin(2*pi*f2*(x+2.5e-9)+p2)*exp(-(x+2.5e-9)/t2)+D;','options',OPS,'problem',{'q','f1','f2'},'coefficients',{'A','k','B','p1','t1','C','p2','t2','D'});
    
    [f0 gof]=fit(fixed_short(strt:end,1),fixed_short(strt:end,2),TYPE,'problem',{q,peak(1),peak(2)});
   

    diffusivity=f0.k;

else
    f0=fit(fixed_short(:,1),fixed_short(:,2),mod_str);
    diffusivity=0;
end

if plotty
    %plot fit to see how well it matches
    figure()
    plot(fixed_short(strt:end,1),f0(fixed_short(strt:end,1)),'b',fixed_short(strt:end,1),fixed_short(strt:end,2),'r');
end

%this is the thermal decay filtering of the signal recorded vs time. This will clean up the DC end of the power spectrum
flat=[fixed_short(strt:end,1) fixed_short(strt:end,2)-f0(fixed_short(strt:end,1))];

if plotty
    figure()
    plot(flat(:,1),flat(:,2),'b-')
end

% Time step info necessary for differentiation and flat padding
tstep=flat(end,1)-flat(end-1,1);

% If option selected, take transform of derivative of recorded signal
% to filter out DC even more than just the background subtraction
if derivative
    d_flat=diff(flat(:,2))/tstep;
    flat=[flat(1:length(d_flat),1) d_flat];
end

% Find the stuff we need to take the spectral profile
num=length(flat(:,1));
fs=num/(flat(end,1)-flat(1,1));
p=18; %magnitude of zero padding to increase resolution in power spectrum
pdsize=2^p-num-2; %more padding = smoother transform

%Only pad on the positive end
pad_val=mean(flat(end-50:end,2));
pad=zeros(pdsize,1);
pad(1:end)=pad_val;
tpad=flat(end,1):tstep:flat(end,1)+(pdsize-1)*tstep;

flat_pad=[flat(:,1) flat(:,2);tpad' pad];

nfft=length(flat_pad(:,2));

%Find the Power Spectral density

%Use a hamming window and a Welchs method. Hamming does the best of the
%ones I've tried and Welch does slightly better than the normal
%periodogram.
[psd freq]=periodogram(flat_pad(:,2),rectwin(nfft),nfft,fs); %periodogram method

%Don't save out DC spike in FFT
amp=sqrt(psd(5:end));
fft=[freq(5:end) amp];

if saveout
    dlmwrite('dat_spec.txt',out);
end

if plotfft
    figure()
    hold on
    plot(freq(5:end),amp,'r');
    xlim([0 1.7e9]);
end

end