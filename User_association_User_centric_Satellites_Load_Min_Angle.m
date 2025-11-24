function [Alpha_association, Satellites_Load] = User_association_User_centric_Satellites_Load_Min_Angle(sr,N_UC,N_max_UTs_per_satellite,Num_Total_satellites,index_satellites,theta_angles_UT,phi_angles_UT,min_Angle_SAT)

Q=size(sr,1);
Alpha_association=cell(Q,1);
Satellites_Load=zeros(Num_Total_satellites,1);
UTs_Associated_to_Satellites=cell(Num_Total_satellites,1);
for qq=1:Q
    N_UC_qq=N_UC;
    sr_qq=sr{qq,1};
    P_qq=length(sr_qq);
    Alpha_association{qq}=zeros(P_qq,1);
    
    if N_UC>P_qq

        N_UC_qq=P_qq;
        
        disp(['Not enough satellites in visibility for the user ', num2str(qq),', the algorithm used ',num2str(N_UC_qq),' serving satellites']);
        
    end
    
    [sr_qq_ordered,indexes_sr]=sort(sr_qq,'ascend');
    
    ind_nn=0;
    
    Satellites_serving_qq=[];
    
    cond_while=1;
    
    while cond_while
        
        ind_nn=ind_nn+1;
        
        candidate_sat=index_satellites{qq,1}(indexes_sr(ind_nn));
        
        cond_angle=1;
        
        elevation_qq=theta_angles_UT{qq}(indexes_sr(ind_nn));

        azimuth_qq=phi_angles_UT{qq}(indexes_sr(ind_nn));
        
        if ~isempty(UTs_Associated_to_Satellites{candidate_sat})
        
            for qq_prime_idx=1:length(UTs_Associated_to_Satellites{candidate_sat})
                
                qq_prime=UTs_Associated_to_Satellites{candidate_sat}(qq_prime_idx);
                
                idx_candidate_sat_qq_prime=find(index_satellites{qq_prime}==candidate_sat);
                
                elevation_qq_prime=theta_angles_UT{qq_prime}(idx_candidate_sat_qq_prime);

                azimuth_qq_prime=phi_angles_UT{qq_prime}(idx_candidate_sat_qq_prime);
                
                gamma_angle=acos( sin(elevation_qq) * sin(elevation_qq_prime) + cos(elevation_qq) * cos(elevation_qq_prime) * cos(azimuth_qq-azimuth_qq_prime));
                
                if (abs(gamma_angle)<min_Angle_SAT)
                    
                    cond_angle=0;
                    
                end
                
            end
            
        end
        
        
        if (Satellites_Load(candidate_sat)<N_max_UTs_per_satellite && cond_angle)
            
            Satellites_serving_qq=[Satellites_serving_qq; indexes_sr(ind_nn)];
            
            Satellites_Load(candidate_sat)=Satellites_Load(candidate_sat)+1;
            
            UTs_Associated_to_Satellites{candidate_sat}=[UTs_Associated_to_Satellites{candidate_sat}; qq];
            
        end
        
        cond_while=(length(Satellites_serving_qq)<N_UC_qq && ind_nn<length(indexes_sr));
        
    end
    
    
    
    Alpha_association{qq,1}(Satellites_serving_qq)=1;
        
  
end
end
