 % find the max segmentation of a biological sequence using CG and MJD
 % method
  function position_genome_max=CGMJD(test_seq)
    score=[];
    score=ntentropy(test_seq);
    negative=find(score<0);
    positive=find(score>0);
    index_position=zeros(1,length(test_seq));
    index_position(negative)=-1;
    index_position(positive)=1;
    % use value to index the change position
    index_change=index_position(1:length(test_seq)-1).*index_position(2:length(test_seq));
    index_position=find(index_change==-1);
    total_number=length(index_position);
    % the position with maxmum MJSD value
    score_position=zeros(1,length(test_seq));
    for zz=1:total_number
        score_position(1,index_position(zz))=mjd(test_seq,index_position(zz));
    end
    [tem position_genome_max]=max(score_position);
 end