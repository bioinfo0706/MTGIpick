 
function label_temp=MSA(Seq,score_nucletide,cg_content,L,sig2,sele,isbound,upstream,downstream)

 %% Path to fast convolution
   if ispc
       addpath([pwd '\Win_CONVNFFT_Folder']);
       %rmpath([pwd '\Unix_CONVNFFT_Folder']);
   elseif isunix
       addpath([pwd '/Unix_CONVNFFT_Folder']);
       %rmpath([pwd '/Win_CONVNFFT_Folder']);
   end
 N=length(Seq);
 % Segment the genome based on the use the score of each window
 [KM,SegmentEnd] = MSS(score_nucletide',L);
 % Enrichment computation
 Enrichment2 = SignificantFoldChange(SegmentEnd,score_nucletide',[],cg_content'); %V2 versus V1 without mappability map
 %% Identify segments based on the scores Enrichment2
 for numb=1:L
    % Detect GIs in each scale
    fprintf(['Detect GIs in ' num2str(numb) 'th scale\n' ]);
    f=[];x=[];
    [f,x]=ksdensity(Enrichment2{numb},'npoints',512);
    CI=FindCI(x,f,sig2);
    index=find(Enrichment2{numb}>CI{1}(2));
    index_01=zeros(1,size(Enrichment2{numb},1));
    index_01(index)=1;
    % sele the begin and end of predicted label and detele some GI with length less than 7
    temp_predict_ttest=index_01;
    % find the star and end position window with value 1
    tem_window_isindex_com=index';
    island_position_window_ttest=zeros(1,2);
    tt=1;
    segment_label=1;
       while segment_label<=length(tem_window_isindex_com) 
              extend_window=0;
              while segment_label+extend_window<length(tem_window_isindex_com)&tem_window_isindex_com(segment_label+1+extend_window)-tem_window_isindex_com(segment_label+extend_window)<2
                    extend_window=extend_window+1;
              end
              island_position_window_ttest(tt,1)=tem_window_isindex_com(segment_label);
              island_position_window_ttest(tt,2)=tem_window_isindex_com(segment_label+extend_window);
              segment_label=segment_label+extend_window+1;
              tt=tt+1;
       end
       for zz=1:tt-1
           if island_position_window_ttest(zz,1)==1
              index_01(island_position_window_ttest(zz,1):island_position_window_ttest(zz,2))=0;
           elseif island_position_window_ttest(zz,2)==size(Enrichment2{numb},1)
              index_01(island_position_window_ttest(zz,1):island_position_window_ttest(zz,2))=0;
%                 elseif island_position_window_ttest(zz,2)-island_position_window_ttest(zz,1)<7
%                     index_01(island_position_window_ttest(zz,1):island_position_window_ttest(zz,2))=0;
           end
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Label the prediction of GIs
   label_predict=zeros(1,length(Seq));
   label_predict(1,1:SegmentEnd{numb}(1)-1)=index_01(1);
   for ii=2:length(index_01)
        label_predict(1,SegmentEnd{numb}(ii-1):SegmentEnd{numb}(ii)-1)=index_01(ii);
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the star and end position window with value 1
   index2=find(label_predict==1);
   tem_window_isindex_com=index2;
   island_position_window_ttest2=zeros(1,2);
   tt=1;
   segment_label=1;
   while segment_label<=length(tem_window_isindex_com) 
         extend_window=0;
          while segment_label+extend_window<length(tem_window_isindex_com)&tem_window_isindex_com(segment_label+1+extend_window)-tem_window_isindex_com(segment_label+extend_window)<2
                extend_window=extend_window+1;
          end
          island_position_window_ttest2(tt,1)=tem_window_isindex_com(segment_label);
          island_position_window_ttest2(tt,2)=tem_window_isindex_com(segment_label+extend_window);
          segment_label=segment_label+extend_window+1;
          tt=tt+1;
   end
   for zz=1:tt-1
       if  island_position_window_ttest2(zz,2)-island_position_window_ttest2(zz,1)<sele
           label_predict(island_position_window_ttest2(zz,1):island_position_window_ttest2(zz,2))=0;
       end
   end         
   %%%%%%%%%%%%%%%%%%%%%%
   label_sequence=label_predict;
   index=find(label_sequence==1);
   tem_window_isindex_com=index;
   island_position_window=zeros(1,2);
   tt=1;
   segment_label=1;
   while segment_label<=length(tem_window_isindex_com) 
          extend_window=0;
          while segment_label+extend_window<length(tem_window_isindex_com)&tem_window_isindex_com(segment_label+1+extend_window)-tem_window_isindex_com(segment_label+extend_window)<2
               extend_window=extend_window+1;
          end
          island_position_window(tt,1)=tem_window_isindex_com(segment_label);
          island_position_window(tt,2)=tem_window_isindex_com(segment_label+extend_window);
          segment_label=segment_label+extend_window+1;
          tt=tt+1;
     end   
    GI_sequence_predict=island_position_window;
    if GI_sequence_predict(1,1)+GI_sequence_predict(1,2)==0
        GI_sequence_predict=[];
    end   
   %% use MJD to detect the boundary
   % find the change position using the ratio between CG and AT
   % the first position detection using CG content and Markov JS divergence
   % Parameters 
   if  isbound==1
       num_GI_predict=size(GI_sequence_predict,1);
       GI_boundary=zeros(num_GI_predict,2);
       for order_GI_predict=1:num_GI_predict
           left_start=max(GI_sequence_predict(order_GI_predict,1)-upstream,1);
           left_end= min(GI_sequence_predict(order_GI_predict,2),GI_sequence_predict(order_GI_predict,1)+ceil(upstream/2));
           right_start=max(GI_sequence_predict(order_GI_predict,1),GI_sequence_predict(order_GI_predict,2)-ceil(downstream/2));
           right_end=min(GI_sequence_predict(order_GI_predict,2)+downstream,length(Seq));    
           if  left_end>length(Seq)
               left_end=length(Seq);
           end 
           if  left_end<1
               left_end=1;
           end 
           test_seq1=Seq(left_start:left_end);
           if right_start<1
              right_start=1;
           end
           if right_start > right_end
              right_length=right_start;
              right_start=right_end;
              right_end=right_length;
           end

          if right_end <= N;
             test_seq2=Seq(right_start:right_end);
             % the boundary
             GI_boundary(order_GI_predict,1)=CGMJD(test_seq1)+left_start;
             GI_boundary(order_GI_predict,2)=CGMJD(test_seq2)+right_start;
             %fprintf(['GI predict: ' num2str(order_GI_predict) '\n']);
          else
             break;
          end
      end
   else
     GI_boundary=GI_sequence_predict;   
 end
%%%%%%%%%%%%%%%%%%%%%%
     label_temp{numb}=GI_boundary; 
 end
end