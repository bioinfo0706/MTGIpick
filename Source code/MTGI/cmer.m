
function [Index Nlettercount signa_str]=cmer(Seq,window,slidelen,k)
% Calculate k-mer of the input data
% user's setting. 
%   Seq-> Input genome
%   window->Size of sliding window
%   slidelen-> Step of sliding windows.
%   k-> Length of k-mer 

%   MTGIpick enables robust identification of genomic islands from a single genome
%   Qi Dai, 20 Apri 2016 


signa_str=cell(1,4^k);
no=1;
N=length(Seq);
N_word=4^k;
T_region=1;
N_region= ceil((length(Seq)-window+1)/slidelen);
merk_c=zeros(N_region,N_word);
NU='ATCG';

% tranlate DNA sequence into matrix data
section=1:slidelen:N-window+1;
Nlettercount=zeros(1,length(section));
Data_t=[];
for i=1:length(section)-1
    subsequence=[];subsequence1=[];
    subsequence=Seq(section(i):section(i)+window);
    subsequence1=subsequence;
    subsequence1(subsequence1=='N')='';
    extend_word=1; 
    subsequence=Seq(section(i):section(i)+window+extend_word-2);
    Nlettercount(1,i)=length(subsequence); 
    Data_t=[Data_t;subsequence];
end

%k-mer
no=1;
if k==2
   for k1=1:4
       for k2=1:4
           a=[];
           a=[NU(k1) NU(k2)];
           signa_str{no}=char(a);
           no=no+1;
       end
   end
elseif k==3
      for k1=1:4
          for k2=1:4
              for k3=1:4
                  a=[];
                  a=[NU(k1) NU(k2) NU(k3)];
                  signa_str{no}=char(a);
                  no=no+1;
              end
          end
       end
elseif k==4   
      for k1=1:4
          for k2=1:4
              for k3=1:4
                     for k4=1:4
                         a=[];
                         a=[NU(k1) NU(k2) NU(k3) NU(k4)];
                         signa_str{no}=char(a);
                         no=no+1;
                     end
                end
           end
       end
elseif k==5
       for k1=1:4
          for k2=1:4
              for k3=1:4
                     for k4=1:4
                         for k5=1:4
                             a=[];
                             a=[NU(k1) NU(k2) NU(k3) NU(k4) NU(k5)];
                             signa_str{no}=char(a);
                             no=no+1;
                         end
                     end
                end
           end
       end
end

Index=zeros(length(section),no-1);

for i=1:length(section)-1
    for j=1:no-1
        Index(i,j)=(length(strfind(Data_t(i,:),signa_str{j}))+1)/(window+4^k);
    end
end

% the last window
    t=1;
    subsequence=[];
    subsequence=Seq(section(length(section)):length(Seq));
    Nlettercount(1,length(section))=length(subsequence);
    if k==2
        for k1=1:4
            for k2=1:4
                a=[];b=[];
                a=[NU(k1) NU(k2)];
                b=findstr(subsequence,a);
                Index(length(section),t)=(length(b)+1)/(window+4^k);
                t=t+1;
            end
        end
    elseif k==3
        for k1=1:4
            for k2=1:4
                for k3=1:4
                    a=[];b=[];
                    a=[NU(k1) NU(k2) NU(k3)];
                    b=findstr(subsequence,a);
                    Index(length(section),t)=(length(b)+1)/(window+4^k);
                    t=t+1;
                end
            end
        end
    elseif k==4
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         a=[];b=[];
                         a=[NU(k1) NU(k2) NU(k3) NU(k4)];
                         b=findstr(subsequence,a);
                         Index(length(section),t)=(length(b)+1)/(window+4^k);
                         t=t+1;
                     end
                end
            end
        end
    elseif k==5
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         for k5=1:4
                             a=[];b=[];
                             a=[NU(k1) NU(k2) NU(k3) NU(k4) NU(k5)];
                             b=findstr(subsequence,a);
                             Index(length(section),t)=(length(b)+1)/(window+4^k);
                             t=t+1;
                         end
                     end
                end
            end
        end
end 
  Nlettercount(1,length(section))=length(subsequence);

end

       