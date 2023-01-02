% Version 1.0: (01/02/2023)
% written by Y. Park

% Yongsung Park & Peter Gerstoft
% MPL/SIO/UCSD
% yongsungpark@ucsd.edu / gerstoft@ucsd.edu
% noiselab.ucsd.edu

% Citation
% Y. Park, P. Gerstoft, and J. H. Lee, "Difference-Frequency MUSIC for DOAs," IEEE Signal Process. Lett. 29, 2612–2616 (2022).

%%
clear; clc;
close all;

% addpath([cd,'/_common'])
% addpath(['../_common'])
errCut = 10; % Maximum RMSE cut-off.

rngNumber = 10000; rng(rngNumber);
dbstop if error;

run_bpdn = 0;

%%
% Environment parameters
% c = 1500;       % speed of sound
c = 343;       % speed of sound

dfreq = 10000;                % difference frequency [Hz]

% ULA-horizontal array configuration
d = .5*(c/dfreq);            % intersensor spacing (based on DF)
% d = 3.75;                  % intersensor spacing [m]
% f = 200; % Hz -> 0.5 lambda = 0.5 * 1500 / 200 = 3.75 [m]

Nsensor = 20;               % number of sensors
q = (0:1:(Nsensor-1))';     % sensor numbering
xq = (q-(Nsensor-1)/2)*d;   % sensor locations

% 5-7.5 times to Difference Frequency (DF)
fmul1 = 5;
fmul2 = 7.5;
ftmp = linspace(fmul1*dfreq,fmul2*dfreq,50);
f = [ftmp,ftmp+dfreq];      % frequency [Hz]
Nfreq = numel(f);

lambda = c./f;   % wavelength


% signal generation parameters
SNR = 20;

% total number of snapshots
Nsnapshot = 50;

% range of angle space
thetalim = [-90 90];

theta_separation = .005;

% Angular search grid
theta = (thetalim(1):theta_separation:thetalim(2))';
Ntheta = length(theta);

% Generate received signal
anglesTrue = [-15; -20];
fprintf(['True DOAs                   :',...
    repmat([' %.4f '],1,numel(anglesTrue)),'\n'],anglesTrue.')

anglesTracks = repmat(anglesTrue,[1,Nsnapshot]);
sinAnglesTracks = sind(anglesTracks); 
Nsource = numel(anglesTrue);

receivedSignalMultiFreq = zeros(Nsensor,Nsnapshot,Nfreq);
Xlist = [];
rPhase = (exp(1i*2*pi*rand(Nsource,Nsnapshot)));
for nfreq = 1:Nfreq
receivedSignal = zeros(Nsensor,Nsnapshot);
e = zeros(Nsensor,Nsnapshot);
for snapshot = 1:Nsnapshot
    source_amp(:,snapshot) = 10*ones(size(anglesTrue));

    Xsource = source_amp(:,snapshot).*exp(1i*2*pi*rand(Nsource,1));    % random phase
    Xlist = [Xlist,Xsource];
   
    % Represenation matrix (steering matrix)
    transmitMatrix = exp( -1i*2*pi/lambda(nfreq)*xq*sinAnglesTracks(:,snapshot).' );
    
    % Received signal without noise
    receivedSignal(:,snapshot) = sum(transmitMatrix*diag(Xsource),2);
    
    % add noise to the signals
    rnl = 10^(-SNR/20)*norm(Xsource);
    nwhite = complex(randn(Nsensor,1),randn(Nsensor,1))/sqrt(2*Nsensor);
    e(:,snapshot) = nwhite * rnl;	% error vector
    receivedSignal(:,snapshot) = receivedSignal(:,snapshot) + e(:,snapshot);

end
receivedSignalMultiFreq(:,:,nfreq) = receivedSignal;
end

% Design/steering matrix (Sensing matrix)
sin_theta = sind(theta);

sensingMatrixFreq = zeros(Nsensor,Ntheta,Nfreq);
for nfreq = 1:Nfreq
    sensingMatrixFreq(:,:,nfreq) = exp(-1i*2*pi/lambda(nfreq)*xq*sin_theta.')/sqrt(Nsensor);
end

% Sensing matrix for DF
sensingMatrixDF = exp(-1i*2*pi/ (c/dfreq) *xq*sin_theta.')/sqrt(Nsensor);


%% Auto-product
APset = zeros(Nsensor,Nsnapshot,Nfreq/2);
for nfreq = 1:Nfreq/2
    receivedSignalTmp(:,:)  = receivedSignalMultiFreq(:,:,nfreq);
    receivedSignalTmp2(:,:) = receivedSignalMultiFreq(:,:,nfreq+Nfreq/2);
    APTmp = receivedSignalTmp2 .* conj(receivedSignalTmp);
    
    APset(:,:,nfreq) = APTmp;
end

%% CBF (ambiguity surface) conventional
PcbfOld = pagemtimes(sensingMatrixFreq,'ctranspose',receivedSignalMultiFreq,'none');
PcbfOld = mean(power(abs(PcbfOld),2),2);
PcbfOld = reshape(PcbfOld,Ntheta,Nfreq);

PcbfOld1 = PcbfOld(:,1:Nfreq/2);
figure, set(gcf,'position',[1,527,560,420])
imagesc(f(1:Nfreq/2)/1e3,theta,10*log10(PcbfOld1./max(PcbfOld1,[],1)))

caxis([-40 0])
set(gca,'fontsize',24,'TickLabelInterpreter','latex', 'YTick',-80:40:80, 'YDir','normal')
set(gca,'position',[0.130857194589887,0.163047626747349,0.774142805410113*2.84/3,0.761952373252651*2.84/3]);
% title('Conventional CBF','interpreter','latex');
xlabel('Lower frequency [kHz]','interpreter','latex');
ylabel('DOA~[$^\circ$]','interpreter','latex');

PcbfOld2 = PcbfOld(:,Nfreq/2+1:Nfreq);
figure, set(gcf,'position',[1,527,560,420])
imagesc(f(Nfreq/2+1:Nfreq)/1e3,theta,10*log10(PcbfOld2./max(PcbfOld2,[],1)))

caxis([-40 0])
set(gca,'fontsize',24,'TickLabelInterpreter','latex','YTick',-80:40:80, 'YDir','normal')
set(gca,'position',[0.130857194589887,0.163047626747349,0.774142805410113*2.84/3,0.761952373252651*2.84/3]);
% title('Conventional CBF','interpreter','latex');
xlabel('Upper frequency [kHz]','interpreter','latex');
% ylabel('DOA~[$^\circ$]','interpreter','latex');
set(gca,'YTickLabel',' ');

%% DF-CBF (ambiguity surface) freq.
PcbfDFfreq = pagemtimes(sensingMatrixDF,'ctranspose',APset,'none');
PcbfDFfreq = mean(power(abs(PcbfDFfreq),2),2);
PcbfDFfreq = reshape(PcbfDFfreq,Ntheta,Nfreq/2);
figure, set(gcf,'position',[560,527,560,420])
imagesc(f(1:Nfreq/2)/1e3,...
    theta,10*log10(PcbfDFfreq./max(PcbfDFfreq,[],1)))

caxis([-40 0])
set(gca,'fontsize',24,'TickLabelInterpreter','latex','YTick',-80:40:80, 'YDir','normal')
% title('DF CBF','interpreter','latex');
set(gca,'position',[0.130857194589887,0.163047626747349,0.774142805410113,0.761952373252651]);
xlabel('Lower frequency [kHz]','interpreter','latex');
ylabel('DOA~[$^\circ$]','interpreter','latex');
% set(gca,'XTickLabel',' ');

%% DF-MUSIC freq.
Neig = Nsource*Nsource;
if Nsnapshot > Nsource*Nsource
PmusicDFfreq = zeros(Ntheta,Nfreq/2);
for nfreq = 1:Nfreq/2
    APsetTmp = reshape(APset(:,:,nfreq),Nsensor,Nsnapshot);
    RzzTmp = APsetTmp*APsetTmp'/Nsnapshot;
    [Rv,Rd] = eig(RzzTmp);
    Rvn = Rv(:, 1:end-Neig);
    PmusicTmp = zeros(numel(theta),1);
    for ii=1:length(theta)
        PmusicTmp(ii) = 1./(sensingMatrixDF(:,ii)'*(Rvn*Rvn')*sensingMatrixDF(:,ii));
    end
    PmusicDFfreq(:,nfreq) = PmusicTmp;
end
PmusicDFfreq = abs(PmusicDFfreq);
figure, set(gcf,'position',[560*2,527,560,420])
imagesc(f(1:Nfreq/2)/1e3,...
    theta,10*log10(PmusicDFfreq./max(PmusicDFfreq,[],1)))

caxis([-40 0])
set(gca,'fontsize',24,'TickLabelInterpreter','latex','YTick',-80:40:80, 'YDir','normal')
% title('DF MUSIC','interpreter','latex');
xlabel('Lower frequency [kHz]','interpreter','latex');
ylabel('DOA~[$^\circ$]','interpreter','latex');
end



%% Auto-Product aligned along columns as they have the same freq. difference
APsetR = reshape(APset,Nsensor,Nsnapshot*Nfreq/2);
RzzTmp = APsetR*APsetR'/(Nsnapshot*Nfreq/2);

%% DF-CBF ALL
PcbfDFall = zeros(numel(theta),1);
for ii=1:length(theta)
    PcbfDFall(ii) = (sensingMatrixDF(:,ii)'*(RzzTmp)*sensingMatrixDF(:,ii));
end
PcbfDFall = abs(PcbfDFall);

[~, Ilocs] = findpeaks(abs(PcbfDFall),'SORTSTR','descend','Npeaks', Nsource);
DoA_error = errorDOAcutoff(theta(Ilocs),anglesTrue,errCut);
disp(['RMSE DF-CBF                 : ',num2str(sqrt(mean(power(DoA_error,2))))])

%% DF-MUSIC ALL
Neig = Nsource;
% RzzTmp = APsetR*APsetR'/(Nsnapshot*Nfreq/2);
[Rv,Rd] = eig(RzzTmp);
Rvn = Rv(:, 1:end-Neig);
PmusicDFall = zeros(numel(theta),1);
for ii=1:length(theta)
    PmusicDFall(ii) = 1./(sensingMatrixDF(:,ii)'*(Rvn*Rvn')*sensingMatrixDF(:,ii));
end
PmusicDFall = abs(PmusicDFall);

[~, Ilocs] = findpeaks(mean(PmusicDFfreq,2),'SORTSTR','descend','Npeaks', Nsource);
DoA_error = errorDOAcutoff(theta(Ilocs),anglesTrue,errCut);
disp(['RMSE Frequency-DF-MUSIC     : ',num2str(sqrt(mean(power(DoA_error,2))))])

[~, Ilocs] = findpeaks(abs(PmusicDFall),'SORTSTR','descend','Npeaks', Nsource);
DoA_error = errorDOAcutoff(theta(Ilocs),anglesTrue,errCut);
disp(['RMSE Time-Frequency-DF-MUSIC: ',num2str(sqrt(mean(power(DoA_error,2))))])

% spRmusic = rmusic_1d(RzzTmp, Neig, 2*pi*d/(c/dfreq));
% DoA_error = errorDOAcutoff(-rad2deg(spRmusic.x_est),anglesTrue,errCut);
% disp(['RMSE Time-Frequency-root-DF-MUSIC: ',num2str(sqrt(mean(power(DoA_error,2))))])



%% DF-MUSIC time, DF-MUSIC-Single snapshot thanks to uniform DF
Neig = Nsource;
if Nfreq/2 > Nsource
PcbfDFsingle = zeros(Ntheta,Nsnapshot);
PmusicDFsingle = zeros(Ntheta,Nsnapshot);
for nsnapshot = 1:Nsnapshot
    APsetTmp = reshape(APset(:,nsnapshot,:),Nsensor,Nfreq/2);
    RzzTmp = APsetTmp*APsetTmp'/(Nfreq/2);
    [Rv,Rd] = eig(RzzTmp);
    Rvn = Rv(:, 1:end-Neig);
    PcbfTmp = zeros(numel(theta),1);
    PmusicTmp = zeros(numel(theta),1);
    for ii=1:length(theta)
        PcbfTmp(ii) = (sensingMatrixDF(:,ii)'*(RzzTmp)*sensingMatrixDF(:,ii));
        PmusicTmp(ii) = 1./(sensingMatrixDF(:,ii)'*(Rvn*Rvn')*sensingMatrixDF(:,ii));
    end
    PcbfDFsingle(:,nsnapshot) = abs(PcbfTmp);
    PmusicDFsingle(:,nsnapshot) = abs(PmusicTmp);
    
    [~, Ilocs] = findpeaks(abs(PcbfTmp),'SORTSTR','descend','Npeaks', Nsource);
    DoA_error_sinCBF(:,nsnapshot) = errorDOAcutoff(theta(Ilocs),anglesTrue,errCut);
    
    [~, Ilocs] = findpeaks(abs(PmusicTmp),'SORTSTR','descend','Npeaks', Nsource);
    DoA_error_sin(:,nsnapshot) = errorDOAcutoff(theta(Ilocs),anglesTrue,errCut);
end
disp(['RMSE Time-DF-CBF   single   : ',num2str(sqrt(mean(mean(power(DoA_error_sinCBF,2)))))])
disp(['RMSE Time-DF-MUSIC single   : ',num2str(sqrt(mean(mean(power(DoA_error_sin,2)))))])

[~, Ilocs] = findpeaks(mean(PmusicDFsingle,2),'SORTSTR','descend','Npeaks', Nsource);
DoA_error = errorDOAcutoff(theta(Ilocs),anglesTrue,errCut);
disp(['RMSE Time-DF-MUSIC          : ',num2str(sqrt(mean(power(DoA_error,2))))])

%% DF-DOA single snapshot contour
PcbfDFsingle = abs(PcbfDFsingle);
figure, set(gcf,'position',[560*1,54,560,420])
imagesc(1:Nsnapshot,theta,10*log10(PcbfDFsingle./max(PcbfDFsingle,[],1)))

caxis([-40 0])
set(gca,'fontsize',24,'TickLabelInterpreter','latex')
title('DF CBF single snapshot','interpreter','latex');
xlabel('Snapshot','interpreter','latex');
ylabel('DOA~[$^\circ$]','interpreter','latex');

PmusicDFsingle = abs(PmusicDFsingle);
figure, set(gcf,'position',[560*2,54,560,420])
imagesc(1:Nsnapshot,theta,10*log10(PmusicDFsingle./max(PmusicDFsingle,[],1)))

caxis([-40 0])
set(gca,'fontsize',24,'TickLabelInterpreter','latex','YTick',-80:40:80, 'YDir','normal')
xlabel('Snapshot','interpreter','latex');
ylabel('DOA~[$^\circ$]','interpreter','latex');
end

%%
figure, set(gcf,'position',[1,54,560,420])
hold on;
plot((PcbfDFall/max(PcbfDFall)),theta,'k','linewidth',1.8)
hold off;

axis([0 1 -90 90]); box on;
set(gca,'fontsize',24,'TickLabelInterpreter','latex','YTick',-80:40:80, 'YDir','normal')
set(gca,'position',[0.130857194589887,0.163047626747349,0.774142805410113,0.761952373252651]);
% title('DF CBF','interpreter','latex');
% ylabel('DOA~[$^\circ$]','interpreter','latex');
xlabel('$P$ [re max]','interpreter','latex');
legend( ['DF-CBF'],...
    'interpreter','latex','location','northeast')
set(gca,'YTickLabel',' ');
% set(gca,'XTickLabel',' ');

%%
figure, set(gcf,'position',[560,54,560,420])
hold on;
plot((mean(PmusicDFfreq,2)/max(mean(PmusicDFfreq,2))),theta,'linewidth',2)
plot((mean(PmusicDFsingle,2)/max(mean(PmusicDFsingle,2))),theta,'g--','linewidth',1.5)
plot((PmusicDFall/max(PmusicDFall)),theta,'r-.','linewidth',1.8)
hold off;

axis([0 1 -90 90]); box on;
set(gca,'fontsize',24,'TickLabelInterpreter','latex','YTick',-80:40:80, 'YDir','normal')
% title('DF MUSIC','interpreter','latex');
% ylabel('DOA~[$^\circ$]','interpreter','latex');
xlabel('$P$ [re max]','interpreter','latex');
legend( ['DF-MUSIC (21)'],...
        ['DF-MUSIC (22)'],...
        ['DF-MUSIC (23)'],...
    'interpreter','latex','location','northeast')
set(gca,'YTickLabel',' ');

%%
% rmpath([cd,'/_common'])
% rmpath(['../_common'])