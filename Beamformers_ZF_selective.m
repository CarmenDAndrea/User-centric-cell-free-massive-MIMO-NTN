function [W] = Beamformers_ZF_selective(index_UT_interest,theta_angles_id, phi_angles_id,theta_angles_eff, phi_angles_eff,...
    Antenna_pos,Gamma_dBi,lambda,Alpha_matrix,index_satellites)
% ZF to only the users that each satellite serves
Q=size(theta_angles_eff,1);
W=cell(Q,1);
N_antenna=size(Antenna_pos,1);
for qq=1:Q
    P_qq=length(theta_angles_eff{qq,1});
    W{qq,1}=zeros(N_antenna,P_qq);
    for pp=1:P_qq
        if Alpha_matrix{qq}(pp)==1
            if qq==index_UT_interest
                theta_angle_pp_qq=theta_angles_id(pp);
                phi_angle_pp_qq=phi_angles_id(pp);
            else
                theta_angle_pp_qq=theta_angles_eff{qq}(pp);
                phi_angle_pp_qq=phi_angles_eff{qq}(pp);
            end
            W_pp_qq=Array_response_generic(Antenna_pos, Gamma_dBi, theta_angle_pp_qq, phi_angle_pp_qq, lambda);
            
            I_int=[];
            sat_index_pp_qq=index_satellites{qq,1}(pp);
            for qq_prime=1:Q
                if qq_prime~=qq
                    sat_index_pp_qq_prime=find(index_satellites{qq_prime}==sat_index_pp_qq);
                    if ~isempty(sat_index_pp_qq_prime)
                        if Alpha_matrix{qq_prime}(sat_index_pp_qq_prime)==1 %selectivity condition: only if the satellites also serves user qq_prime
                            if qq_prime==index_UT_interest
                                theta_angle_pp_qq_prime=theta_angles_id(sat_index_pp_qq_prime);
                                phi_angle_pp_qq_prime=phi_angles_id(sat_index_pp_qq_prime);
                            else
                                theta_angle_pp_qq_prime=theta_angles_eff{qq_prime}(sat_index_pp_qq_prime);
                                phi_angle_pp_qq_prime=phi_angles_eff{qq_prime}(sat_index_pp_qq_prime);
                            end
                            W_pp_qq_prime=Array_response_generic(Antenna_pos, Gamma_dBi, theta_angle_pp_qq_prime, phi_angle_pp_qq_prime, lambda);
                            
                            I_int=[I_int, W_pp_qq_prime];
                        end
                    end
                end
            end
            
            if ~isempty(I_int)
                
                I_orth=orth(I_int);
                
                W_pp_qq_new=(eye(N_antenna)-I_orth*I_orth')*W_pp_qq;
                
            else
                W_pp_qq_new=W_pp_qq;
            end
            W{qq}(:,pp)=conj(W_pp_qq_new)/norm(W_pp_qq_new);
            
            
        end
    end
end
end

