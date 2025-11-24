function [paths_energy,Delays,Doppler_frequencies,theta_angles_SAT_id, phi_angles_SAT_id,theta_angles_SAT_eff, phi_angles_SAT_eff,...
    theta_angles_UT_id, phi_angles_UT_id,theta_angles_UT_eff, phi_angles_UT_eff] = Data_Elaboration_p681(sr,channelLMS,del,fs,elSAT,azSAT,elUT,azUT,f,f_star)
%Assumes that the first user is the reference one and its actual position
%is referred to user 2
Q=size(sr,1);
Q_eff=Q-1;
lambda=3e8/f;
paths_energy=cell(Q_eff,1);
Delays=cell(Q_eff,1);
Doppler_frequencies=cell(Q_eff,1);
theta_angles_SAT_eff=cell(Q_eff,1);
phi_angles_SAT_eff=cell(Q_eff,1);
theta_angles_UT_eff=cell(Q_eff,1);
phi_angles_UT_eff=cell(Q_eff,1);

% Ideal position is the one of user 1 (reference user)
elSAT_ref=elSAT{1,1};
azSAT_ref=azSAT{1,1};
elUT_ref=elUT{1,1};
azUT_ref=azUT{1,1};
[theta_angles_SAT_id, phi_angles_SAT_id] = Geographic_to_polarConversion(elSAT_ref, azSAT_ref);
[theta_angles_UT_id, phi_angles_UT_id] = Geographic_to_polarConversion(elUT_ref, azUT_ref);

del_ref=del{1,1};

fs_ref=fs{1,1};

% Remove reference user from data

sr_eff=sr(2:end);

del_eff=del(2:end);

fs_eff=fs(2:end);

elSAT_eff=elSAT(2:end);

azSAT_eff=azSAT(2:end);

elUT_eff=elUT(2:end);

azUT_eff=azUT(2:end);

channelLMS_eff=channelLMS(2:end);


for qq=1:Q_eff
    
    if qq==1
        Delays_noNorm_qq=abs(del_eff{1,1}-del_ref); %relative delay with respect to the reference user
    
        Delays{qq,1}=Delays_noNorm_qq-min(Delays_noNorm_qq);
        
        Doppler_frequencies{qq,1}=fs{1,1}*f/f_star-fs_ref*f/f_star; %Doppler frequencies normalized to the frequencies f_star
    else
        Delays{qq,1}=del_eff{qq,1}; %absolute delay
        
        Doppler_frequencies{qq,1}=fs_eff{qq-1,1}*f/f_star; %absolute Doppler
    end
    sr_qq=sr_eff{qq,1};
    elSAT_qq=elSAT_eff{qq,1};
    azSAT_qq=azSAT_eff{qq,1};
    elUT_qq=elUT_eff{qq,1};
    azUT_qq=azUT_eff{qq,1};
    chann_qq=channelLMS_eff{qq};
    
    paths_energy{qq}=(lambda./(4*pi*sr_qq)).*chann_qq;
    
    [theta_angles_SAT_eff{qq,1}, phi_angles_SAT_eff{qq,1}] = Geographic_to_polarConversion(elSAT_qq, azSAT_qq);

    [theta_angles_UT_eff{qq,1}, phi_angles_UT_eff{qq,1}] = Geographic_to_polarConversion(elUT_qq, azUT_qq);
end

end     