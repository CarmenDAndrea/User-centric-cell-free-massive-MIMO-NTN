function [Psi_matrices] = Psi_matrix_OFDM_reduced_size(index_UT_interest,N,M,Delays, Doppler_frequencies,Delta_f,T)

P_index_UT_interest=size(Delays{index_UT_interest},1);

Psi_matrices=cell(P_index_UT_interest,1);
for pp=1:P_index_UT_interest
    Psi_matrices{pp,1}=zeros(N*M,N*M);
end

I_N=eye(N);

mPrimeVec = (0:M-1)';

for pp = 1:P_index_UT_interest
    
    multCoeff = exp(-1i*2*pi*mPrimeVec*Delta_f*Delays{index_UT_interest}(pp)).*squeeze(sum(exp(1i*2*pi*mPrimeVec.*reshape((Doppler_frequencies{index_UT_interest}(pp)/Delta_f+mPrimeVec-(0:M-1))/M,[1 M M])),1));
    
    for n = 0:N-1
        for m = 0:M-1
            psi_temp=1/M*exp(1i*2*pi*n*T*Doppler_frequencies{index_UT_interest}(pp)).*multCoeff(:,m+1);
            psi_kron=kron(psi_temp,I_N(:,n+1));
            Psi_matrices{pp,1}(m*N+n+1,:)=transpose(psi_kron);
        end
    end
    
end
end

