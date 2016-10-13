function  [ScoringList  ks_position_whole  ks_ascending_whole]=WholeComputeScore(mer_whole,method,feature_set,eye_window)
   
    ScoringList=[];
%% select the main window of whole genome
    core_window_whole=[];
    variance_window_whole=[];
    variance_window_whole=var(mer_whole')';
    densitiesSig1_var_whole=[];
    densitiesSig2_var_whole=[];
    [densitiesSig1_var_whole,densitiesSig2_var_whole,xmaxlist,ymaxlist]=CalculateDensities(variance_window_whole, 1,0);% '5sigma' or '3sigma'    
    ScoringList_var_whole=[];
    ScoringList_var_whole=ComputeScore(densitiesSig1_var_whole,densitiesSig2_var_whole,variance_window_whole);
    core_window_whole=mer_whole;
    core_window_whole(find(ScoringList_var_whole==1),:)=[];

 %% select features using kertosis and skewness
    num_feature_whole=size(mer_whole,2);
    ks_value_whole=zeros(2,num_feature_whole);
    for i1=1:num_feature_whole
         ks_value_whole(1,i1)= kurtosis(mer_whole(:,i1)); % 
         ks_value_whole(2,i1) = skewness(mer_whole(:,i1));  % 
    end
         ks_ascending_whole=zeros(2,num_feature_whole);
         ks_position_whole=zeros(2,num_feature_whole);
    for i1=1:2
        [ks_ascending_whole(i1,:) ks_position_whole(i1,:)]=sort(ks_value_whole(i1,:),'descend');
    end   
     Chosen_feature=ks_position_whole(1,1:feature_set);
     
  %% compute score  
    if strmatch(method,'kernal')
       %% compute the score using kernal density function
        densitiesSig1_whole=[];
        densitiesSig2_whole=[];
        [densitiesSig1_whole,densitiesSig2_whole,xmaxlist,ymaxlist]=CalculateDensities(core_window_whole, 3,4);% '5sigma' or '3sigma'
        ScoringList=ComputeScore(densitiesSig1_whole,densitiesSig2_whole,mer_whole,Chosen_feature);%largest
    else
       %% compute the score using eye and main t-test method
        ScoringList=eandmvalue(mer_whole,core_window_whole,eye_window,Chosen_feature);% largest
    end
end

 
    
    
    
    
    
    