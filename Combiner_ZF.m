function [V] = Combiner_ZF(index_UT_interest, theta_angles_id, phi_angles_id,theta_angles_eff, phi_angles_eff,...
    Antenna_pos,Gamma_dBi,lambda,Alpha_matrix)
Q=size(theta_angles_eff,1);
N_antenna=size(Antenna_pos,1);
V=cell(Q,1);
for qq=1:Q
    P_qq=length(theta_angles_eff{qq,1});
    V{qq,1}=zeros(N_antenna,P_qq);
    for pp_star=1:P_qq
        if Alpha_matrix{qq}(pp_star)==1
            if qq==index_UT_interest
                theta_angle_pp_star_qq=theta_angles_id(pp_star);
                phi_angle_pp_star_qq=phi_angles_id(pp_star);
            else
                theta_angle_pp_star_qq=theta_angles_eff{qq}(pp_star);
                phi_angle_pp_star_qq=phi_angles_eff{qq}(pp_star);
            end
            a_p_star_qq=Array_response_generic(Antenna_pos, Gamma_dBi, theta_angle_pp_star_qq, phi_angle_pp_star_qq, lambda);
            
            I_int=[];
            
            for pp_star2=1:P_qq
                if Alpha_matrix{qq}(pp_star2)==1
                    if pp_star2~=pp_star
                        if qq==index_UT_interest
                            theta_angle_pp_star2_qq=theta_angles_id(pp_star2);
                            phi_angle_pp_star2_qq=phi_angles_id(pp_star2);
                        else
                            theta_angle_pp_star2_qq=theta_angles_eff{qq}(pp_star2);
                            phi_angle_pp_star2_qq=phi_angles_eff{qq}(pp_star2);
                        end
                        a_UT_pp_qq=Array_response_generic(Antenna_pos, Gamma_dBi, theta_angle_pp_star2_qq, phi_angle_pp_star2_qq, lambda);
                        
                        I_int=[I_int, a_UT_pp_qq];
                    end
                end
                
            end
            
            if ~isempty(I_int)
                
                I_orth=orth(I_int);
                
                V_pp_star_qq=(eye(N_antenna)-I_orth*I_orth')*a_p_star_qq;
                
            else
                V_pp_star_qq=a_p_star_qq;
            end

              V{qq}(:,pp_star)=conj(V_pp_star_qq)/norm(V_pp_star_qq); % no compensation of the channel gain
            
        end
    end
end

end


