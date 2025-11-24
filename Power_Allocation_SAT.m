function [Powers_DL] = Power_Allocation_SAT(Pt_SAT_dBW,Alpha_matrix,M,N)

Pt_SAT=10^(Pt_SAT_dBW/10);

P=size(Alpha_matrix,1);

Q=size(Alpha_matrix,2);

Powers_DL=zeros(P,Q);

for pp=1:P
    
    if sum(Alpha_matrix(pp,:))~=0
    
        Powers_DL(pp,:)=Pt_SAT.*Alpha_matrix(pp,:)/M/N; 
        % uniformly distributed power between the users served by the pp-th satellite and the symbols
    
    end
end

end

