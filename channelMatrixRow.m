function [PsiICI,PsiISI,lengthISI] = channelMatrixRow(tau,nu,N,M,T,DeltaF,Dir2,calcDir2)
kVec = 0:N-1;
lPrimeVec = 0:M-1;
lVec = transpose(0:M-1);
if(calcDir2)
    Dir2 = sin(pi*(-lVec+lPrimeVec+tau*M*DeltaF))./sin(pi/M*(-lVec+lPrimeVec+tau*M*DeltaF))...
        .*exp(1i*pi*(-lVec+lPrimeVec+tau*M*DeltaF)*(M-1)/M);
    Dir2(isnan(Dir2)) = M;
end
lICI = 0:M-1-ceil(tau/(T/M));
lISI = M-1-floor(tau/(T/M)):M-1;
expTermICI = exp(1i*2*pi*nu*(lICI/(M*DeltaF)));
expTermISI = exp(1i*2*pi*nu*(lISI/(M*DeltaF)-T));
if (length(expTermISI)+length(expTermICI) == M+1)
    expTerm = [expTermICI(1:end) expTermISI(2:end)];
    lengthISI = length(expTermISI(2:end));
else  
    expTerm = [expTermICI expTermISI];
    lengthISI = length(expTermISI);
end
if max(size(Dir2))~=max(size(expTerm)) %Modifica per risolvere il problema di size mismatch in qualche scenario
    diffDim=max(size(Dir2))-max(size(expTerm));
    expTerm=[expTerm, ones(1,  diffDim)];
end
Dir2Exp = 1/(N*M)*Dir2.*expTerm;
Dir2ExpR = reshape(Dir2Exp,[M*M 1]);
expK =exp(-1i*2*pi*kVec/N);
PsiISI = Dir2ExpR(M*M-lengthISI*M+1:end)*expK;
PsiICI = repmat(Dir2ExpR(1:M*M-lengthISI*M),1,N);
end