function [Powers_DL] = Uniform_power_Allocation(Pt_SAT_dBW,Alpha_association,M,N,index_satellites,Num_Total_satellites)


%% Alpha_matrix generation

Q=size(index_satellites,1);

Alpha_matrix_full=zeros(Num_Total_satellites,Q);

for qq=1:Q   
    Alpha_qq=Alpha_association{qq,1};
    P_qq=length(Alpha_qq);
    for pp=1:P_qq
        if Alpha_qq(pp)==1
            Alpha_matrix_full(index_satellites{qq,1}(pp),qq)=1;
        end
    end 
end

%% Power distribution
Powers_DL_full=zeros(Num_Total_satellites,Q);

Pt_SAT=10^(Pt_SAT_dBW/10);

for sat=1:Num_Total_satellites
    if sum(Alpha_matrix_full(sat,:))~=0
        Powers_DL_full(sat,:)=Pt_SAT/sum(Alpha_matrix_full(sat,:)).*Alpha_matrix_full(sat,:)/M/N; 
        % uniformly distributed power between the users served by the pp-th satellite and the symbols
    end
end

%% Save only the required data

Powers_DL=cell(Q,1);

for qq=1:Q
    Powers_DL{qq,1}=Powers_DL_full(index_satellites{qq,1},qq);
end

end

