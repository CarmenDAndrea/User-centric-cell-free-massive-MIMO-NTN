function [Psi_matrices] = Psi_matrix_OTFS_reduced_size(index_UT_interest,N,M,Delays, Doppler_frequencies,Delta_f,T,definePulses,pulseTx,pulseRx)

P_index_UT_interest=size(Delays{index_UT_interest},1);

Psi_matrices=cell(P_index_UT_interest,1);
for pp=1:P_index_UT_interest
    Psi_matrices{pp,1}=zeros(N*M,N*M);
end

for pp = 1:P_index_UT_interest
    
    PsiDec = complex(zeros(M*M,N));
    [PsiICI,PsiISI,lengthISI] = ...
        channelMatrixRow(Delays{index_UT_interest}(pp),Doppler_frequencies{index_UT_interest}(pp),N,M,T,Delta_f,[],1);               % The real channel matrix, without approximations for the interference carriers
    transition_ICI_ISI = M*M-lengthISI*M;
    % Determine the index of the part of the matrix which belongs to ISI or ISI computation
    PsiDec(1:transition_ICI_ISI,:)     = PsiICI;                                                       % Save the part of Psi relative to ICI
    PsiDec(transition_ICI_ISI+1:end,:) = PsiISI;                                                       % Save the part of Psi relative ti ISI
    Dir1Shift = Dir1ShiftCalc(Doppler_frequencies{index_UT_interest}(pp),T,N,M);
    Dir1ShiftR = repelem(Dir1Shift,1,M);
    if(definePulses == 0)                                                                              % If no pulses are specified, use the default rectangular pulses
        Psi_matrices{pp,1}=repmat(reshape(PsiDec,[M M*N]),N,1).*transpose(Dir1ShiftR);
    else
        Psi_matrices{pp,1}=channelMatrixPulses(Psi_matrices{pp,qq},Delays{index_UT_interest}(pp),Doppler_frequencies{index_UT_interest}(pp),...
            N,M,T,Delta_f,pulseTx,pulseRx);                                           % Use the specified pulses with their own channel matrix
    end
    
end

end

