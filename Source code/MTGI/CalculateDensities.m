function  [densitiesSig1,densitiesSig2,xmaxlist,ymaxlist]=CalculateDensities(data, filtration,k)

densitiesSig1=cell(1,4^k);
densitiesSig2=cell(1,4^k);
xmaxlist=[];
ymaxlist=[];
%outlayers=ClusterFiltration(data,filtration);
outlayers=ExtremeFiltration(data,filtration);
for l=1:4^k
    densities1=cell(1,length(data));
    densities2=cell(1,length(data));
	xmax=0;
	ymax=0;
% 	for k =1:length(data)
	    if length(outlayers{l})>0
                tee= data(:,l);
                te=outlayers{l};
                tee(te)=[];
				[f,x]=ksdensity(tee,'npoints',512);				
        else
                [f,x]=ksdensity(data(:,l),'npoints',512);
        end
		%	names(densities)[k]=paste("chromosome",k,sep="")
	    if (max(x)>xmax)
              xmax=max(x);
        end
	    if (max(f)>ymax)
              ymax=max(f);
        end
%         densities1{k}=x; 
%         densities2{k}=f; 
%     end
	xmaxlist(l)=xmax;
	ymaxlist(l)=ymax;
	densitiesSig1{l}=x;
    densitiesSig2{l}=f;
  end
end
