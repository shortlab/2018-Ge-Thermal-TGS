function [fft] = TGS_phase_fft(pos_file,neg_file,grat,psd_out)
%   For this use only, generate filtered fft for input data
%   pos_file: positive phase data pos_file
%   neg_file: negative phase data neg_file
%   grat: grating, either calibrated or estimated, in um
%   psd: if psd=1, save out power spectrum else, save out fft magnitude

if nargin<4
    psd_out=0;
end

%Boolean options for plotting and processing
derivative=0;
plotty=0;
plotfft=0;
saveout=0;

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

fixed_short(:,1)=fixed_short(:,1)-fixed_short(1,1); %slide peak back to 0 to fit more easily

%Take a fit to a purely erfc() decay to subtract a decay profile before
%taking fft
q=2*pi/(grat*10^(-6));
LB=[0 0 0];
UB=[10 5*10^-4 0.1];
ST=[0.05 5*10^-5 0];
OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
TYPE=fittype('A.*erfc(q*sqrt(k*(x)))+c;','options',OPS,'problem','q','coefficients',{'A','k','c'});

[f0,~]=fit(fixed_short(:,1),fixed_short(:,2),TYPE,'problem',q);

if plotty
    %plot fit to see how well it matches
    figure()
    plot(fixed_short(:,1),f0(fixed_short(:,1)),'b',fixed_short(:,1),fixed_short(:,2),'r');
end

%this is the thermal decay filtering of the signal recorded vs time. This will clean up the DC end of the power spectrum
flat=[fixed_short(:,1) fixed_short(:,2)-f0(fixed_short(:,1))];

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
[psd,freq]=periodogram(flat_pad(:,2),rectwin(nfft),nfft,fs); %periodogram method

%Don't save out DC spike in FFT/PSD
if psd_out
    amp=psd(5:end);
else
    amp=sqrt(psd(5:end));
end

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