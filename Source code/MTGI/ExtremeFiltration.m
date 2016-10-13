function outlayers=ExtremeFiltration(data,filtration)
%outlayers=cell(1,length(data));
%if (strcmp(filtration,none)){}
%if strcmp(filtration,'5sigma')
%	for (chromosome=1 :length(data))
        [nrow,ncol]=size(data);
	    chromOutlayers=cell(1,ncol);
		for tetranucleotide=1:ncol
		     TetraOutlayers1=(find(data(:,tetranucleotide)>(mean(data(:,tetranucleotide)+filtration*std(data(:,tetranucleotide))))))';
		     TetraOutlayers2=(find(data(:,tetranucleotide)<(mean(data(:,tetranucleotide)-filtration*std(data(:,tetranucleotide))))))';
		     chromOutlayers{tetranucleotide}=[TetraOutlayers1,TetraOutlayers2];
        end
	    outlayers=chromOutlayers;
end
