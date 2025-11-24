%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the Channel Matrix for Arbitrary Pulses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Psi = channelMatrixPulses(Psi,tau,nu,N,M,T,DeltaF,pulseTx,pulseRx)
kVec = 0:N-1;
kPrimeVec = transpose(0:N-1);
lPrimeVec = 0:M-1;
lVec = transpose(0:M-1);
Dir2 = sin(pi*(-lVec+lPrimeVec+tau*M*DeltaF))./sin(pi/M*(-lVec+lPrimeVec+tau*M*DeltaF))...
    .*exp(1i*pi*(-lVec+lPrimeVec+tau*M*DeltaF)*(M-1)/M);                                                           % Delay Dirichlet kernel function
Dir2(isnan(Dir2)) = M;
Dir1 = sin(pi*(nu*N*T+kPrimeVec-kVec))./sin(pi/N*(nu*N*T+kPrimeVec-kVec))...
    .*exp(1i*pi*(nu*N*T+kPrimeVec-kVec)*(N-1)/N);                                                                  % Doppler Dirichlet kernel function
Dir1(isnan(Dir1)) = N;
Dir1(abs(Dir1)<1e-5) = 0;
gTx = transpose(pulseTx(lPrimeVec/M))/(N*M);                                                                       % Sampled transmitted pulse with normalization
if(sum(isnan(gTx)))
    NanIdx = isnan(gTx);
    gTx(NanIdx) = pulseTx((lPrimeVec(NanIdx)+0.000001)/M)/(N*M);
end
lPrimeExp = transpose(exp(1i*2*pi*nu*lPrimeVec/(M*DeltaF)));
nInt = 10;                                                                                                         % Number of interfering pulses
interfVec = -nInt:nInt;                                                                                            % Take into account the difference (n-n') which counts the number of pulses interfering in each sampling point
interfVec = reshape(interfVec,[1 1 length(interfVec)]);
notNanPulse = conj(pulseRx((lPrimeVec*T/M-interfVec*T+tau)/T));
if(sum(sum(isnan(notNanPulse))))
    notNanPulse(isnan(notNanPulse)) = conj(pulseRx((lPrimeVec(NanIdx)*T/M-0*T+tau+0.000001)/T));
end
notNanPulse(notNanPulse == -Inf) = 0;
curlyTerm =  exp(-1i*2*pi*(kPrimeVec+nu*N*T).*interfVec/N).*notNanPulse;                                           % Term taking into account the interference from all pulses in a certain sampling instant. It represents the summation over n and n'
curlyTerm = sum(curlyTerm,3);                                                                                      % Sum over all possible differences (n-n')
multFact = gTx.*lPrimeExp.*Dir2;
multFact(abs(multFact)<1e-5) = 0;
for kPrime = 0:N-1
    multFact1 = multFact.*curlyTerm(kPrime+1,:);
    Psi(kPrime*M+(0:M-1)+1,:) = reshape(reshape(reshape(multFact1,[M*M 1])*Dir1(kPrime+1,:),[M M N]),[M M*N]);
end
end