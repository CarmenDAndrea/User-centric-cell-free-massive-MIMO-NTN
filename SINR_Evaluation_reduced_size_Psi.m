function [SINR] = SINR_Evaluation_reduced_size_Psi(index_UT_interest,Alpha_matrix,Delta_matrices,Psi_matrices, A_matrix_SAT, A_matrix_UT, W, V, Rho_val, D_matrices,N_matrices,Powers, noise_variance,index_satellites)
NM=size(Psi_matrices{1},1);
D_qq=D_matrices;
N_qq=N_matrices;
Q=size(Alpha_matrix,1);
P_index_UT_interest=size(Psi_matrices,1);
SINR=zeros(NM,1);

S_qq=find(Alpha_matrix{index_UT_interest});
D_qq_H=D_matrices';
for ll=1:NM
    Sum_num_qq_ll=0;
    
    Sum_Int_qq_ll_1=zeros(NM,1);
    
    Sum_Int_qq_ll_2=zeros(NM,1);
    
    Sum_Int_qq_ll_3=zeros(NM,Q);
    
    for pp_star_idx=1:length(S_qq)
        pp_star=S_qq(pp_star_idx);
        mu_p_star_p_star_q_q=Rho_val{index_UT_interest}(pp_star)*V{index_UT_interest}(:,pp_star).'*A_matrix_UT{index_UT_interest}(:,pp_star)*A_matrix_SAT{index_UT_interest}(:,pp_star).'*W{index_UT_interest}(:,pp_star);
        Sum_num_qq_ll=Sum_num_qq_ll+Delta_matrices{pp_star}(ll,ll)*sqrt(Powers{index_UT_interest}(pp_star))*mu_p_star_p_star_q_q*D_qq_H(ll,:)*Psi_matrices{pp_star}(:,ll);
        for ii=1:NM
            if ii~=ll
                Sum_Int_qq_ll_1(ii,1)=Sum_Int_qq_ll_1(ii,1)+Delta_matrices{pp_star}(ll,ll)*sqrt(Powers{index_UT_interest}(pp_star))*mu_p_star_p_star_q_q*D_qq_H(ll,:)*Psi_matrices{pp_star}(:,ii);
            end
            
            for pp=1:P_index_UT_interest
                if pp~=pp_star
                    if (Alpha_matrix{index_UT_interest}(pp)==1)
                        mu_p_star_p_q_q=Rho_val{index_UT_interest}(pp)*V{index_UT_interest}(:,pp_star).'*A_matrix_UT{index_UT_interest}(:,pp)*A_matrix_SAT{index_UT_interest}(:,pp).'*W{index_UT_interest}(:,pp);
                        Sum_Int_qq_ll_2(ii,1)=Sum_Int_qq_ll_2(ii,1)+Delta_matrices{pp_star}(ll,ll)*sqrt(Powers{index_UT_interest}(pp))*mu_p_star_p_q_q*D_qq_H(ll,:)*Psi_matrices{pp}(:,ii);
                    end
                end
                
                sat_index_pp_qq=index_satellites{index_UT_interest}(pp);
                
                for qq_prime=1:Q
                    if qq_prime~=index_UT_interest
                        sat_index_pp_qq_prime=find(index_satellites{qq_prime}==sat_index_pp_qq);
                        if ~isempty(sat_index_pp_qq_prime)
                            if Alpha_matrix{qq_prime}(sat_index_pp_qq_prime)==1
                                mu_p_star_p_q_q_prime=Rho_val{index_UT_interest}(pp)*V{index_UT_interest}(:,pp_star).'*A_matrix_UT{index_UT_interest}(:,pp)*A_matrix_SAT{index_UT_interest}(:,pp).'*W{qq_prime}(:,sat_index_pp_qq_prime);
                                Sum_Int_qq_ll_3(ii,qq_prime)=Sum_Int_qq_ll_3(ii,qq_prime)+Delta_matrices{pp_star}(ll,ll)*sqrt(Powers{qq_prime}(sat_index_pp_qq_prime))*mu_p_star_p_q_q_prime*D_qq_H(ll,:)*Psi_matrices{pp}(:,ii);
                            end
                        end
                    end
                    
                    
                end
                
            end
        end
        
    end
    
    SINR_num_qq_ll=abs(Sum_num_qq_ll)^2;
    
    SINR_den_qq_ll=sum(abs(Sum_Int_qq_ll_1).^2+abs(Sum_Int_qq_ll_2).^2)+ ...
        sum(sum(abs(Sum_Int_qq_ll_3).^2))+noise_variance*norm(D_qq_H(ll,:)*N_matrices)^2;
    
    if SINR_den_qq_ll==0
        SINR(ll)=0;
    else
        SINR(ll)=SINR_num_qq_ll/SINR_den_qq_ll;
    end
    
    
end
SINR=real(SINR);

end

