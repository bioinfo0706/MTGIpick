
function  [window_isindex_whole_kernal island_position_window]=GIwindowposition(Tem_score,significant)
   
     %kernal density function
     f=[];x=[];
     [f,x]=ksdensity(Tem_score,'npoints',512);	
     CI=FindCI(x,f,significant);%
     window_isindex_whole_kernal=find(Tem_score>CI{1}(2));
     tem_window_isindex_com=window_isindex_whole_kernal;
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
end


  