function [Alpha_association] = User_association_User_centric_distances(sr,N_UC)

Q=size(sr,1);
Alpha_association=cell(Q,1);

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
    
    Satellites_serving=indexes_sr(1:N_UC_qq);
    
    Alpha_association{qq,1}(Satellites_serving)=1;
        
  
end

end
