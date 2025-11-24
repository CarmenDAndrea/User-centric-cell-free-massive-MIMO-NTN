function [index_satellites,sr,channLMS,del,fs,elSAT,azSAT,elUT,azUT] = Remove_NaN_from_Data_p681(sr_NaN,channLMS_NaN,del_NaN,fs_NaN,elSAT_NaN,azSAT_NaN,elUT_NaN,azUT_NaN)
Q=size(sr_NaN,2);

index_satellites=cell(Q,1);
sr=cell(Q,1);
del=cell(Q,1);
fs=cell(Q,1);
elSAT=cell(Q,1);
azSAT=cell(Q,1);
elUT=cell(Q,1);
azUT=cell(Q,1);
channLMS=cell(Q,1);

numChIt=size(channLMS_NaN,3);

for qq=1:Q
    
    sr_NaN_qq=sr_NaN(:,qq);
    index_satellites{qq,1}=find(~isnan(sr_NaN_qq));
    channLMS{qq,1}=zeros(length(index_satellites{qq,1}),numChIt);

    if qq==2
        length1=length(index_satellites{1});
        length2=length(index_satellites{2});
        min_length=min(length1,length2);
       if length2~=length1
           index_satellites{2}=index_satellites{1};
          
       elseif sum(index_satellites{2}(1:min_length)~=index_satellites{1}(1:min_length))>0 
           index_satellites{2}=index_satellites{1};
          
       end
    end
    
    sr{qq,1}=sr_NaN(index_satellites{qq,1},qq);

    channLMS{qq,1}=squeeze(channLMS_NaN(index_satellites{qq,1},qq,:));
    
    del{qq,1}=del_NaN(index_satellites{qq,1},qq);
    
    fs{qq,1}=fs_NaN(index_satellites{qq,1},qq);
    
    elSAT{qq,1}=elSAT_NaN(index_satellites{qq,1},qq);
    
    azSAT{qq,1}=azSAT_NaN(index_satellites{qq,1},qq);
    
    elUT{qq,1}=elUT_NaN(index_satellites{qq,1},qq);
    
    azUT{qq,1}=azUT_NaN(index_satellites{qq,1},qq);

end

end