function   ScoringList=ComputeScore(densityWMM,densitiesSig2,data,ChosenTetra) %densityWMM=densitiesSig1

[nrow,ncol]=size(densityWMM);
%     ScoringList=cell(1,length(data));
% 	for k =1:length(data)
%        ScoringList{k}=zeros(size(data{k},1),1);
%     end
ScoringList=zeros(size(data,1),1);
tempVec=zeros(size(data,1),ncol); 
	for l=1:ncol
% 		for k =1:length(data)
		    CI=FindCI(densityWMM{l},densitiesSig2{l},0.05);%[0.05, 0.025, 0.01]
			for AlphaNumber=1:length(CI)
               sz1=find(data(:,l)<CI{AlphaNumber}(1));
               sz2=find(data(:,l)>CI{AlphaNumber}(2));
               ssz=zeros(size(data,1),1);
               ssz(sz1)=[1];
               ssz(sz2)=[1];
% 			   tempVec(:,l)=tempVec(:,l)+ssz.*data(:,l);
 			   tempVec(:,l)=tempVec(:,l)+ssz;
            end
%		    ScoringList=[ScoringList;tempVec];	
%         end	
    end
if nargin>3
   ScoringList=sum(tempVec(:,ChosenTetra)')';
else
   if ncol==1
     ScoringList=tempVec;
   else
    ScoringList= sum(tempVec')';
   end
   
end
end
