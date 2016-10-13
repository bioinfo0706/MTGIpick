
function Enrichment = SignificantFoldChange(SegmentEnd,V1,UM,V2,pth)

%% Scenarios
BackgroundFlag = 0;
VersusFlag = 0;
if exist('UM','var');
    if ~isempty(UM);
        BackgroundFlag = 1;
    end
end
if exist('V2','var');
    if ~isempty(V2);
        VersusFlag = 1;
        %versus is in beta
    end
end

%% Parameters
minLsize1 = 0;
minLsize2 = 10;
th = 1e-5;
if ~exist('pth','var');
    pth = 1e-6;
end
maxEFC = 10;                        %only in versus mode (1=fold_change)
foldchange_or_additive = 1;         %only in versus mode (1=fold_change)

%% Loop scales depending on scenario
if VersusFlag==0&BackgroundFlag==0;
    
    MLV = max(V1);
    V1 = V1./MLV;
    UM = ones(size(V1));
    
    CL = cumsum(UM);
    CV = cumsum(V1);
    p = CV(end)/CL(end);
    
    L = length(SegmentEnd);
    Enrichment = cell(L,1);
    
    for nn = 1:L;
%         fprintf(1,'.');
        
       %% Creating variables
        SL = length(SegmentEnd{nn});
        Segments = cat(1,[1 SegmentEnd{nn}(1:SL-1)+1],SegmentEnd{nn});
        Le   = (CL(Segments(2,:))-CL(Segments(1,:))+UM(Segments(1,:)));
        O    = (CV(Segments(2,:))-CV(Segments(1,:))+V1(Segments(1,:)));
        G2 = (max(Le,minLsize1));
        G2r = (max(Le,minLsize2));
        G1 = CV(end)*ones(size(O));
        T = CL(end)*ones(size(O));
        
        %% Computing enrichment
%         Enrichment{nn} = EFC(T,G1,G2,O,pth,th)';
        Enrichment{nn} = EFC_inflatedvariance(T,G1,G2,G2r,O,pth,th);
        
    end
elseif VersusFlag==0&BackgroundFlag==1;
    
    %normalize signals to one
    MLV = max(UM);
    UM = UM./MLV;
    MLV = max(V1);
    V1 = V1./MLV;
    %make sure that the signal is never larger than the background
    V1 = min(V1,UM);
    
    CL = cumsum(UM);
    CV = cumsum(V1);
    p = CV(end)/CL(end);
    
    L = length(SegmentEnd);
    Enrichment = cell(L,1);
    
    for nn = 1:L;
%         fprintf(1,'.');
        
        %% Creating variables
        SL = length(SegmentEnd{nn});
        Segments = cat(1,[1 SegmentEnd{nn}(1:SL-1)+1],SegmentEnd{nn});
        Le   = (CL(Segments(2,:))-CL(Segments(1,:))+UM(Segments(1,:)));
        O    = (CV(Segments(2,:))-CV(Segments(1,:))+V1(Segments(1,:)));
        G2 = (max(Le,minLsize1));
        G2r = (max(Le,minLsize2));
        G1 = CV(end)*ones(size(O));
        T = CL(end)*ones(size(O));
        
        %% Computing enrichment
%         Enrichment{nn} = EFC(T,G1,G2,O,pth,th)';

          Enrichment{nn} = EFC_inflatedvariance(T,G1,G2,G2r,O,pth,th);
    end
    
elseif VersusFlag==1&BackgroundFlag==0;
    
    MLV = max(V1);
    V1 = V1./MLV;
    UM = ones(size(V1));
    V2 = V2./MLV;
    
    CL = cumsum(UM);
    CV = cumsum(V1);
    CR = cumsum(V2);
    
    L = length(SegmentEnd);
    Enrichment = cell(L,1);
    
    for nn = 1:L;
%         fprintf(1,'.');
        
        %% Creating variables
        SL = length(SegmentEnd{nn});
        Segments = cat(1,[1 SegmentEnd{nn}(1:SL-1)+1],SegmentEnd{nn});
        Le   = (CL(Segments(2,:))-CL(Segments(1,:))+UM(Segments(1,:)));
        O    = (CV(Segments(2,:))-CV(Segments(1,:))+V1(Segments(1,:)));
        Rr    = (CR(Segments(2,:))-CR(Segments(1,:))+V2(Segments(1,:)));
        
        G2 = (max(Le,minLsize1));
        G2r = (max(Le,minLsize2));
        
        %% Computing enrichment
        Enrichment{nn} = EFC_versus_iv(G2,G2r,O,Rr,pth,th,maxEFC,foldchange_or_additive);
    end
    
elseif VersusFlag==1&BackgroundFlag==1;
    
    MLV = max(UM);
    V1 = V1./MLV;
    UM = UM./MLV;
    V2 = V2./MLV;
    CL = cumsum(UM);
    CV = cumsum(V1);
    CR = cumsum(V2);
    
    L = length(SegmentEnd);
    Enrichment = cell(L,1);
    
    for nn = 1:L;
%         fprintf(1,'.');
        
        %% Creating variables
        SL = length(SegmentEnd{nn});
        Segments = cat(1,[1 SegmentEnd{nn}(1:SL-1)+1],SegmentEnd{nn});
        Le   = (CL(Segments(2,:))-CL(Segments(1,:))+UM(Segments(1,:)));
        O    = (CV(Segments(2,:))-CV(Segments(1,:))+V1(Segments(1,:)));
        Rr    = (CR(Segments(2,:))-CR(Segments(1,:))+V2(Segments(1,:)));
        
        G2 = (max(Le,minLsize1));
        G2r = (max(Le,minLsize2));
        
        %% Computing enrichment
        Enrichment{nn} = EFC_versus_iv(G2,G2r,O,Rr,pth,th,maxEFC,foldchange_or_additive);
    end
    
end
end
