function [D_matrices,N_matrices] = Practical_LMMSE_Detection_MultiAntenna_UTs(index_UT_interest,Alpha_matrix, Delta_matrices,Psi_matrices, A_matrix_SAT, A_matrix_UT, W, V, Rho_val,Powers, noise_variance)

NM=size(Psi_matrices{1},1);

D_matrices=zeros(NM,NM);

N_UT=size(V{index_UT_interest}(:,1),1);

N_matrices=zeros(NM,NM*N_UT);

A_qq=zeros(NM,NM);

N_qq=zeros(NM,NM*N_UT);

S_qq=find(Alpha_matrix{index_UT_interest});

for pp_star_idx=1:length(S_qq)
    pp=S_qq(pp_star_idx);
    mu_val=Rho_val{index_UT_interest}(pp)*V{index_UT_interest}(:,pp).'*A_matrix_UT{index_UT_interest}(:,pp)*A_matrix_SAT{index_UT_interest}(:,pp).'*W{index_UT_interest}(:,pp);
    A_qq=A_qq+Delta_matrices{pp}*sqrt(Powers{index_UT_interest}(pp))*mu_val*Psi_matrices{pp};
    N_qq=N_qq+kron(Delta_matrices{pp},V{index_UT_interest}(:,pp).');
end

D_matrix_qq_H=A_qq'*inv(A_qq*A_qq'+noise_variance*(N_qq*N_qq'));

D_matrices(:,:)=D_matrix_qq_H';

N_matrices(:,:)=N_qq;


end

