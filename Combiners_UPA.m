function [V] = Combiners_UPA(index_UT_interest, theta_angles_id, phi_angles_id,theta_angles_eff, phi_angles_eff,...
    Antenna_pos,Gamma_dBi,lambda,Alpha_matrix)
Q=size(theta_angles_eff,1);
N_antenna=size(Antenna_pos,1);
V=cell(Q,1);
for qq=1:Q
    P_qq=length(theta_angles_eff{qq,1});
    V{qq,1}=zeros(N_antenna,P_qq);
    for pp=1:P_qq
        if Alpha_matrix{qq}(pp)==1
            if qq==index_UT_interest
                theta_angle_pp_qq=theta_angles_id(pp);
                phi_angle_pp_qq=phi_angles_id(pp);
            else
                theta_angle_pp_qq=theta_angles_eff{qq}(pp);
                phi_angle_pp_qq=phi_angles_eff{qq}(pp);
            end
            a_pp_qq=Array_response_generic(Antenna_pos, Gamma_dBi, theta_angle_pp_qq, phi_angle_pp_qq, lambda);
            V{qq}(:,pp)=conj(a_pp_qq)/norm(a_pp_qq); % no compensation of the channel gain
        end
    end
end
end

