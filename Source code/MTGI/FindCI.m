function   CI=FindCI(densitiesSig1,densitiesSig2,alpha)
CI=cell(1,length(alpha));
dx=mean(diff(densitiesSig1)); % calculate differences between predicted values
dx=dx/(sum(densitiesSig2)*dx); % repair differences so everything sums prettily to 1
cumDens=(densitiesSig2)*dx; 	%compute cumulative distribution function from density (this and step after this)
cumSum=cumsum(cumDens); % I dont like to compute cumsum every time so I save this rather than call it again (change from previous version)
	for i=1:length(alpha)% cycle for every alpha
	    lowerCI=densitiesSig1(min(find(cumSum>(alpha(i)/2)))); % finds lower CI
	    upperCI=densitiesSig1(min(find(cumSum>(1-alpha(i)/2)))); % finds upper CI
	    CI{i}=[lowerCI,upperCI]; % save to list
%    upperCI=ksdensity(densitiesSig1 ,1-alpha(i)/2 ,'function','icdf'); 
%    lowerCI=ksdensity(densitiesSig1,alpha(i)/2  ,'function','icdf');   
%     CI{i}=[lowerCI,upperCI];
    end
end
