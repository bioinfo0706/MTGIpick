function [score_nucletide ks_position_whole ks_ascending_whole]=ISTLFS(Seq,G_signatures,window,slidelen,sig,ww,interation_whole,feature_size,eye_window,each_window_length)
% To reveal potentially dramatic changes at a smaller scale, 
% we propose an iteration of small-scale t-tests with large-scale feature selection (IST-LFS) 
% to quantify the compositional differences of a region from the host
% user¡¯s setting. 
%   Seq-> Input genome
%   G_signatures->Genomic signatures
%   Sig-> Standard deviation of the mean of the window scores to select windows whose scores are large enough to
%       be considered statistically significant.
%   Interation_whole-> Periods of time that are repeated to select windows whose scores are large enough to 
%       be considered statistically significant.
%   feature_size-> Size of selected features by the proposed kurtosis
%   each_window_length-> Length of the windows
%   MTGIpick enables robust identification of genomic islands from a single genome
%   Qi Dai, 20 Apri 2016 

%     G_signatures=mer{1};
%     slidelen=1000;
%     sig=0.05;
%     ww=1;
%     interation_whole=5;
%     feature_size=256;
%     eye_window=5;
    
    %% select the main window of whole genome  
    sig_len=length(sig);
    label_whole_kernal=zeros(sig_len,interation_whole,size(G_signatures,1));
    label_whole_ttest=zeros(sig_len,interation_whole,size(G_signatures,1));

    % combinate larger windows
    [row column]=size(G_signatures);
    for width_window=ww:ww  
        % compute frequency of larger windows
        temp_mer_width=zeros(row-width_window+1,column);
        if width_window==1
            temp_mer_width=G_signatures;
        else
           for order_feature=1:size(G_signatures,2)
                slide_value=zeros(width_window,width_window+size(G_signatures,1)-1);
                for zz=1:width_window
                    slide_value(zz,:)=[zeros(1,width_window-zz) G_signatures(:,order_feature)' zeros(1,zz-1)];
                end
                slide_value(:,1:(width_window-1))=[];
                slide_value(:,(size(G_signatures,1)-(width_window-1)+1):size(G_signatures,1))=[];
                temp_mer_width(:,order_feature)=mean(slide_value)';
           end
        end
    end
    
    for sig_order=1:sig_len
        fprintf([' Compute whole score with Significant: ' num2str(sig(sig_order)) '\n' ]);
        % interation of process  
        tem_mer_interation_kernal=[];tem_mer_interation_ttest=[];
        tem_mer_interation_kernal=temp_mer_width;
        tem_mer_interation_ttest=temp_mer_width;
        N=length(Seq);
        section=1:slidelen:N-window+1;
        label_window_kernal=1:length(section);
        label_window_ttest=1:length(section);
        GI_predict_Interation=zeros(1,length(section));
        Label_window_whole=zeros(2,size(G_signatures,1));
        interation_score=cell(1,interation_whole);%get score of each interation
        interation_index=cell(1,interation_whole);%get index of each interation place
        interation_index2=cell(1,interation_whole);%get index of original place
        for interation_time=1:interation_whole
            % compute score using kernal density or ttest
            %feature_size=4^k;
            ScoringList_kernal_temp=[];ScoringList_ttest_temp=[];
            [ScoringList_ttest_temp ,ks_position_whole, ks_ascending_whole]=WholeComputeScore(tem_mer_interation_ttest,'eandm',feature_size,eye_window);
            ks_position_whole(2,:)=[];
            ks_ascending_whole(2,:)=[];
            % find the significant windows
            position_window_whole=cell(1,2);
            [position_window_whole{2} island_position_window_whole_ttest]=GIwindowposition(ScoringList_ttest_temp,sig(sig_order));
     
            % label of predicted window
            for zzz=2:2
                temp_len_p=length(position_window_whole{zzz});
                for zzzz=1:temp_len_p
                    Label_window_whole(zzz,label_window_ttest(position_window_whole{zzz}(zzzz):position_window_whole{zzz}(zzzz)+width_window-1))=1;
                end
            end
            label_whole_ttest(sig_order,interation_time,:)=Label_window_whole(2,:);
      
            % update the frequency of windows
            interation_index{interation_time}=position_window_whole{2};
            interation_score{interation_time}=ScoringList_ttest_temp(position_window_whole{2});
            tem_mer_interation_ttest(position_window_whole{2},:)=[];
            label_window_ttest(position_window_whole{2})=[];
            fprintf(['Interation: ' num2str(interation_time) '\n' ]);
      end % interation
 end % significant

    score_index=1:length(section);
    for ii= 1:interation_whole
        interation_index2{ii}=score_index(interation_index{ii});
        score_index(interation_index{ii})=[];
    end
    score_our=zeros(length(section)-width_window+1,1);
    for ii= 1:interation_whole-1
        score_our(interation_index2{ii})=interation_score{ii};
    end
    score_our(find(score_our==0))=ScoringList_ttest_temp;

    %% set label to each bases of sequence
    score_nucletide=[];
    label_sequence=zeros(1,length(Seq));
    tem_position=cumsum(each_window_length);
    num_window_total=length(tem_position);
    label_sequence(1,1:tem_position(1,1))=score_our(1);
    for zz=2:(num_window_total-width_window+1)
        label_sequence(1,tem_position(1,zz-1)+1:tem_position(1,zz))=score_our(zz);
    end
   score_nucletide=label_sequence;  
end