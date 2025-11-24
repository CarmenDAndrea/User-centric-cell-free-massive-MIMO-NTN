function [Alpha_association, Satellites_Load] = User_association_User_centric_distances_Satellites_Load(sr,N_UC,N_max_UTs_per_satellite,Num_Total_satellites,index_satellites)

Q=size(sr,1);
Alpha_association=cell(Q,1);
Satellites_Load=zeros(Num_Total_satellites,1);
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
        
        if Satellites_Load(candidate_sat)<N_max_UTs_per_satellite
            
            Satellites_serving_qq=[Satellites_serving_qq; indexes_sr(ind_nn)];
            
            Satellites_Load(candidate_sat)=Satellites_Load(candidate_sat)+1;
            
        end
        
        cond_while=(length(Satellites_serving_qq)<N_UC_qq && ind_nn<length(indexes_sr));
        
    end
    
    
    
    Alpha_association{qq,1}(Satellites_serving_qq)=1;
        
  
end
end
