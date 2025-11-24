function [Satellites_Load] = Evaluate_Satellites_Load(Alpha_association,index_satellites,Num_Total_satellites)


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


Satellites_Load=sum(Alpha_matrix_full,2);

end

