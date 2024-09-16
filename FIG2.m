% Version 1.0: (01/02/2023)
% written by Yongsung Park

% Yongsung Park & Peter Gerstoft
% MPL/SIO/UCSD
% yongsungpark@ucsd.edu / gerstoft@ucsd.edu
% noiselab.ucsd.edu

% Citation
% Y. Park, P. Gerstoft, and J. H. Lee, “Difference-Frequency MUSIC for DOAs,” IEEE Signal Process. Lett. 29, 2612–2616 (2022).
% https://doi.org/10.1109/LSP.2022.3230365

% Atomic norm minimization (ANM) implementation is also available.
% Y. Park and P. Gerstoft, “Difference Frequency Gridless Sparse Array Processing,” IEEE Open J. Signal Process. 5, 914–925 (2024).
% https://doi.org/10.1109/OJSP.2024.3425284

%%
clear; clc;
close all;

% addpath([cd,'/_common'])
% addpath(['../_common'])
errCut = 10; % Maximum RMSE cut-off.

rngNumber = 1; rng(rngNumber);
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
anglesTrue = [-40; 35];
fprintf(['True DOAs        :',...
    repmat([' %.4f '],1,numel(anglesTrue)),'\n'],anglesTrue.')

anglesTracks = repmat(anglesTrue,[1,Nsnapshot]);
anglesTracks(1,:) = anglesTracks(1,1) - 2*anglesTracks(1,1)./(1+exp(-.1*(-Nsnapshot/2:-Nsnapshot/2+Nsnapshot-1)));
anglesTracks(2,:) = anglesTracks(2,1) - 1.00*(0:Nsnapshot-1)';
anglesTracks(:,46) = anglesTracks(:,47);

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


%% Auto-Product aligned along columns as they have the same freq. difference
APsetR = reshape(APset,Nsensor,Nsnapshot*Nfreq/2);

%% MUSIC-Single snapshot thanks to multi-frequency
Neig = Nsource;
if Nfreq/2 > Nsource
PmusicDFsingle = zeros(Ntheta,Nsnapshot);
for nsnapshot = 1:Nsnapshot
    APsetTmp = reshape(APset(:,nsnapshot,:),Nsensor,Nfreq/2);
    RzzTmp = APsetTmp*APsetTmp'/(Nfreq/2);
    [Rv,Rd] = eig(RzzTmp);
    Rvn = Rv(:, 1:end-Neig);
    PmusicTmp = zeros(numel(theta),1);
    for ii=1:length(theta)
        PmusicTmp(ii) = 1./(sensingMatrixDF(:,ii)'*(Rvn*Rvn')*sensingMatrixDF(:,ii));
    end
    PmusicDFsingle(:,nsnapshot) = PmusicTmp;
    
    [~, Ilocs] = findpeaks(abs(PmusicTmp),'SORTSTR','descend','Npeaks', Nsource);
    DoA_error_sin(:,nsnapshot) = errorDOAcutoff(theta(Ilocs),anglesTracks(:,nsnapshot),errCut);    
end
disp(['RMSE single-MUSIC: ',num2str(sqrt(mean(mean(power(DoA_error_sin,2)))))])

PmusicDFsingle = abs(PmusicDFsingle);
figure,
set(gcf,'position',[180,200,1180,420])
set(gca,'position',[0.157,0.163,0.748,0.762])
imagesc(1:Nsnapshot,theta,10*log10(PmusicDFsingle./max(PmusicDFsingle,[],1)))

caxis([-40 0])
set(gca,'fontsize',24,'TickLabelInterpreter','latex','YTick',-80:40:80, 'YDir','normal')
% title('DF MUSIC single snapshot','interpreter','latex');
xlabel('Snapshot','interpreter','latex');
ylabel('DOA~[$^\circ$]','interpreter','latex');
% set(gca,'YTickLabel',' ');
end


%%
% rmpath([cd,'/_common'])
% rmpath(['../_common'])