function [A_matrix_SAT,A_matrix_UT] = Array_response_Evaluation_SAT_UT(Antenna_pos_SAT, Antenna_pos_UT,theta_angles_SAT_eff,...
    phi_angles_SAT_eff,theta_angles_UT_eff,phi_angles_UT_eff,lambda, Gamma_dBi_SAT, Gamma_dBi_UT)
Q=size(theta_angles_SAT_eff,1);
NS=size(Antenna_pos_SAT,1);
N_UT=size(Antenna_pos_UT,1);

A_matrix_SAT=cell(Q,1);
A_matrix_UT=cell(Q,1);


for qq=1:Q
    P_qq=length(theta_angles_SAT_eff{qq,1});
    A_matrix_SAT{qq,1}=zeros(NS,P_qq);
    A_matrix_UT{qq,1}=zeros(N_UT,P_qq);
    for pp=1:P_qq
        A_matrix_SAT{qq,1}(:,pp)=Array_response_generic(Antenna_pos_SAT, Gamma_dBi_SAT, theta_angles_UT_eff{qq,1}(pp), phi_angles_UT_eff{qq,1}(pp), lambda);
        A_matrix_UT{qq,1}(:,pp)=Array_response_generic(Antenna_pos_UT, Gamma_dBi_UT, theta_angles_SAT_eff{qq,1}(pp), phi_angles_SAT_eff{qq,1}(pp), lambda);
    end
end

end

