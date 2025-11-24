function [channLMS,totDopplerFreq] = p681LMS_Channel(fCarrier,fShift,numChIt,elSAT)

maxMobSpeed = (130)*1000/3600; % km/h -> m/s
minMobSpeed = (50)*1000/3600; % km/h -> m/s

P=size(fShift,1);
Q=size(fShift,2);

cLight=3e8; %m/s

% p681LMS Channel (narrowband)
comPar = struct;
comPar.CarrierFrequency = fCarrier;
comPar.MobileAltitude = 0;
comPar.RandomStream = "Global stream";
comPar.NumSinusoids = 48;
mobileMaxDoppler = maxMobSpeed*comPar.CarrierFrequency/cLight;
comPar.SampleRate = round((max(fShift(:))+mobileMaxDoppler)*11); % the sampling rate must be at least ten times the max Doppler
channLMS = zeros(P,Q,numChIt);
totDopplerFreq = zeros(P,Q);
for q = 1:Q   % users
    for p = 1:P  % satellites
        if ~isnan(elSAT(p,q))
        comPar.ElevationAngle = elSAT(p,q);   % from file
        comPar.MobileSpeed = (maxMobSpeed-minMobSpeed)*rand()+minMobSpeed;  % max speed uniformly distr. in the range (min-max)
        satelliteDopplerShift = fShift(p,q);  % residual satellite Doppler
        % % channel object parameters
        ntnNBChan = p681LMSChannel;
        ntnNBChan.SampleRate = comPar.SampleRate;
        ntnNBChan.CarrierFrequency = comPar.CarrierFrequency;
        ntnNBChan.ElevationAngle = comPar.ElevationAngle;
        ntnNBChan.MobileSpeed = comPar.MobileSpeed;
        ntnNBChan.SatelliteDopplerShift = satelliteDopplerShift;
        ntnNBChan.RandomStream = comPar.RandomStream;
        ntnNBChan.Environment = "Highway";
        ntnNBChan.ChannelFiltering = false;
        ntnNBChan.NumSamples = comPar.SampleRate;
        ntnNBChan.AzimuthOrientation = randi([0 360]);
        ntnNBChan.FadingTechnique = "Sum of sinusoids";
        ntnNBChan.NumSinusoids = comPar.NumSinusoids;
        [channTmp,~,~] = ntnNBChan();
        channLMS(p,q,:) = abs(channTmp(round(linspace(1,comPar.SampleRate,numChIt)))); % 1 sec time-series sampled given the number of iterations per snapshot
        release(ntnNBChan); % reset the channel to generate a new realization for each link
        totDopplerFreq(p,q) = satelliteDopplerShift + cosd(ntnNBChan.AzimuthOrientation)*cosd(comPar.ElevationAngle)*comPar.MobileSpeed*comPar.CarrierFrequency/cLight;
        % p681ChannelInfo = info(ntnNBChan)
        end
    end
end




% COMMENTI
% totDopplerFreq contiene la nuova matrice con i valori di Doppler complessivi per ogni canale utente-satellite da usare per generare Psi
% channLMS è la matrice con i valori di ampiezza del canale (solo modulo perché nella fase qui comparirebbe il Doppler di cui noi teniamo conto a parte) con cui moltiplicare la potenza ricevuta in ogni canale utente-satellite

end

