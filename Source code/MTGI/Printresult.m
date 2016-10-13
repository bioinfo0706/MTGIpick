
function label_temp=Printresult(Seq,file_name,feature_size,ks_position_whole,ks_ascending_whole,signa_str,L,label_temp,xunhuan)

% index1=findstr(infile,'.');
% out=infile(1:index1-1);
if xunhuan==0
   out1=strcat(file_name,'_signa','.txt');
   out2=strcat(file_name,'_Total_SPGIs','.gff3');
   out3=strcat(file_name,'_Each_SPGIs','.txt');
   out4=strcat(file_name,'_html','.html');
else
   out1=strcat(file_name,'_signa','_sequence',num2str(xunhuan),'.txt');
   out2=strcat(file_name,'_Total_SPGIs','_sequence',num2str(xunhuan),'.gff3');    
   out3=strcat(file_name,'_Each_SPGIs','_sequence',num2str(xunhuan),'.txt');
   out4=strcat(file_name,'_html','_sequence',num2str(xunhuan),'.html'); 
end

FD1=fopen(out1,'w');
FD31=fopen(out2,'w');
FD=fopen(out3,'w');
FD4=fopen(out4,'w');

N=length(Seq);
% Print the genomic signatures
%fprintf(FD1,'%s\n',Name);
fprintf(['Print genomic signatures \n' ]);
fprintf(FD1,'NO        kmer        kurtosis score\n');
no=1;
format short g;
for hh=1:(feature_size)
   fprintf(FD1,'%d        %s        %d\n',no,signa_str{ks_position_whole(no)},ks_ascending_whole(no));
   no=no+1;
end
fprintf(FD1,'\n');
 fclose(FD1);

%%%%%%%%%%%%%%%%Label all of the bases of the predicted genome
label_all=cell(1,L);
for xh=1:L
      label_sequence=zeros(1,length(Seq));
      for zz=1:size(label_temp{xh},1)
          label_sequence(1,label_temp{xh}(zz,1):label_temp{xh}(zz,2))=1;
      end
      label_all{xh}= label_sequence;
end

%Sum all the predicted value
label_allseq=label_all{1};
for numb=2:L
    label_allseq=label_allseq+label_all{numb};
end

% find the star and end position window with same conserved score
largeone=max(label_allseq);
 xposition=[];
for overlap=0:largeone
    tem_island=find(label_allseq==overlap)';
    if isempty(tem_island)==0       
       island_position=zeros(1,2);
       tt=1;
       segment_label=1;
       while segment_label<=length(tem_island) 
             extend_window=0;
             while segment_label+extend_window<length(tem_island)&tem_island(segment_label+1+extend_window)-tem_island(segment_label+extend_window)<2
                   extend_window=extend_window+1;
             end
             island_position(tt,1)=tem_island(segment_label);
             island_position(tt,2)=tem_island(segment_label+extend_window);
             segment_label=segment_label+extend_window+1;
             tt=tt+1;
       end
       num_island=size(island_position,1);
       px=zeros(num_island,3);
       for i=1:num_island
           px(i,:)=[island_position(i,1) island_position(i,2) overlap];
       end
       xposition{overlap+1}=px;
    end
end

% calculate the conserved scores of each predicted GIs
% caulcuate the predict GI in all scales
canshu=[];
for ijk=1:largeone+1
    canshu=[canshu;xposition{1,ijk}];
end
[canshu1,pos] = sort(canshu(:,1)); 
canshu = canshu(pos,:);
for ijk=1:size(canshu,1)
    canshu(ijk,4)=canshu(ijk,2)-canshu(ijk,1)+1;
end
%detailed caulcuate the predict GI areas of each scales
allcanshu=[];
if canshu(1,3)==0
    i=1;
    for ijk=1:length(xposition{1,1}(:,1))-1
        allcanshu(i,:)=xposition{1,1}(ijk,:);
        xa=[];xb=[];xc=[];temp=[];
        xa=xposition{1,1}(ijk,2)+1;
        xb=xposition{1,1}(ijk+1,1)-1;
        xc=1;
        temp=[xa,xb,xc];
        allcanshu(i+1,:)=temp; 
        i=i+2;
    end
    allcanshu(i,:)=xposition{1,1}(length(xposition{1,1}(:,1)),:);
    if allcanshu(i,2) < N
        xa=[];xb=[];xc=[];temp=[];
        xa=xposition{1,1}(ijk+1,2)+1;
        xb=N;xc=1;
        temp=[xa,xb,xc];
        allcanshu(i+1,:)=temp;
    end
else
    xa=[];xb=[];xc=[];temp=[];
    xa=1;
    xb=xposition{1,1}(1,1)-1;
    xc=1;
    temp=[xa,xb,xc];
    allcanshu(1,:)=temp;
    i=2;
    for ijk=1:length(xposition{1,1}(:,1))-1
        allcanshu(i,:)=xposition{1,1}(ijk,:);
        xa=[];xb=[];xc=[];temp=[];
        xa=xposition{1,1}(ijk,2)+1;
        xb=xposition{1,1}(ijk+1,1)-1;
        xc=1;
        temp=[xa,xb,xc];
        allcanshu(i+1,:)=temp; 
        i=i+2;
    end
    allcanshu(i,:)=xposition{1,1}(length(xposition{1,1}(:,1)),:);
    if allcanshu(i,2) < N
        xa=[];xb=[];xc=[];temp=[];
        xa=xposition{1,1}(ijk+1,2)+1;
        xb=N;xc=1;
        temp=[xa,xb,xc];
        allcanshu(i+1,:)=temp;
    end
end

for ijk=1:size(allcanshu,1)
    allcanshu(ijk,4)=allcanshu(ijk,2)-allcanshu(ijk,1)+1;
end    
%calculate the cg_content 
for cg = 1:1:length(allcanshu(:,1))
    if allcanshu(cg,3) > 0
        start_seq=allcanshu(cg,1);
        end_seq=allcanshu(cg,2);
        seq_str=Seq(start_seq:end_seq);
        ab=[];
        ab=findstr(seq_str,'C');
        c_content=length(ab);
        ab=[];
        ab=findstr(seq_str,'G');
        g_content=length(ab);
        allcanshu(cg,5)=(c_content+g_content)/length(seq_str);
    else
        allcanshu(cg,5)=0;
    end
end

for cg = 1:1:length(canshu(:,1))
    if canshu(cg,3) > 0
        start_seq=canshu(cg,1);
        end_seq=canshu(cg,2);
        seq_str=Seq(start_seq:end_seq);
        ab=[];
        ab=findstr(seq_str,'C');
        c_content=length(ab);
        ab=[];
        ab=findstr(seq_str,'G');
        g_content=length(ab);
        canshu(cg,5)=(c_content+g_content)/length(seq_str);
    else
        canshu(cg,5)=0;
    end
end
%print to the 'a b c d' parmeters
a=allcanshu(:,4);
a_index=allcanshu(:,3);
c=allcanshu(:,5);
b=cell(length(a),4);
d=cell(length(a),1);
mm=1;ijk1=1;
weizhi=find(canshu(:,3)==0);
if canshu(1,3)>0
    temp=canshu(1:weizhi(1)-1,:);
    b{mm,1}=temp(:,4)';
    b{mm,2}=temp(:,3)';
    atemp=temp(:,1:2)';
    btemp=atemp(:);
    b{mm,4}=btemp';
    d{mm,1}=temp(:,5)';
    mm=mm+1;ijk1=1;
    for ijk=2:size(allcanshu,1)-1
        if allcanshu(ijk,3)>0
            temp=canshu(weizhi(ijk1)+1:weizhi(ijk1+1)-1,:);
            b{mm,1}=temp(:,4)';
            b{mm,2}=temp(:,3)';
            atemp=temp(:,1:2)';
            btemp=atemp(:);
            b{mm,4}=btemp';
            d{mm,1}=temp(:,5)';
            ijk1=ijk1+1;
        end
        mm=mm+1;
    end
    if canshu(size(canshu,1),3) > 0
        temp=canshu(weizhi(ijk1)+1:length(canshu),:);
        b{mm,1}=temp(:,4)';
        b{mm,2}=temp(:,3)';
        atemp=temp(:,1:2)';
        btemp=atemp(:);
        b{mm,4}=btemp';
        d{mm,1}=temp(:,5)';
    end
else
    for ijk=1:size(allcanshu,1)-1
        if allcanshu(ijk,3)>0
            temp=canshu(weizhi(ijk1)+1:weizhi(ijk1+1)-1,:);
            b{mm,1}=temp(:,4)';
            b{mm,2}=temp(:,3)';
            atemp=temp(:,1:2)';
            btemp=atemp(:);
            b{mm,4}=btemp';
            d{mm,1}=temp(:,5)';
            ijk1=ijk1+1;
        end
        mm=mm+1;
    end
    if canshu(size(canshu,1),3) > 0
        temp=canshu(weizhi(ijk1)+1:length(canshu),:);
        b{mm,1}=temp(:,4)';
        b{mm,2}=temp(:,3)';
        atemp=temp(:,1:2)';
        btemp=atemp(:);
        b{mm,4}=btemp';
        d{mm,1}=temp(:,5)';
    end
end

%%%%%%%%%%%%%%%%%print out31 and print out3%%%%%%%%%%%%%%%%%%%%%
str1 = '##gff-version 3';
str2 = '##sequence-region';
str3 = 'Scale';
str4 = 'MTGIpick';
str5 = 'genomic_island';
str6 = '+';
str7 = '.';
str8 = 'ID=GI';
%print the out3
for num_m=1:length(label_all)
    fprintf(FD,'%s %d\n',str3,num_m);
    tem_island=find(label_all{num_m}==1);
    island_position=zeros(1,2);
    tt=1;
    segment_label=1;
    while segment_label<=length(tem_island) 
          extend_window=0;
          while segment_label+extend_window<length(tem_island)&tem_island(segment_label+1+extend_window)-tem_island(segment_label+extend_window)<2
                extend_window=extend_window+1;
          end
          island_position(tt,1)=tem_island(segment_label);
          island_position(tt,2)=tem_island(segment_label+extend_window);
          segment_label=segment_label+extend_window+1;
          tt=tt+1;
    end
    for is=1:length(island_position(:,1))
        island_position(is,3)=island_position(is,2)-island_position(is,1)+1;
    end
    fprintf(FD,'%s\n',str1);
    fprintf(FD,'%s %s %d %d\n',str2,file_name,1,length(Seq));
    gi=1;
    for is=1:length(island_position(:,1))
        each_start=island_position(is,1);
        each_end=island_position(is,2);
        each_length=island_position(is,3);
        if each_start > 0 && each_end > 0
            ID=strcat(str8,num2str(gi));
            fprintf(FD,'%s	%s	%s	%d	%d	%s	%s	%s	%s\n',file_name,str4,str5,each_start,each_end,'1',str6,str7,ID);
            gi=gi+1;
        else
            fprintf(FD,'\n');
        end
    end
    fprintf(FD,'\n');
end
%print the out31

if sum(label_allseq) > 0
    fprintf(FD31,'%s\n',str1);
    fprintf(FD31,'%s %s %d %d\n',str2,file_name,1,length(Seq));
    gi=1;
    position=find(canshu(:,3)>0);
    for ijk=1:1:length(position)-1
        ID=strcat(str8,num2str(gi));
        fprintf(FD31,'%s	%s	%s	%d	%d	%d	%s	%s	%s\n',file_name,str4,str5,canshu(position(ijk),1),canshu(position(ijk),2),canshu(position(ijk),3),str6,str7,ID);
        gi=gi+1;
    end
    ID=strcat(str8,num2str(gi));
    fprintf(FD31,'%s	%s	%s	%d	%d	%d	%s	%s	%s',file_name,str4,str5,canshu(position(length(position)),1),canshu(position(length(position)),2),canshu(position(length(position)),3),str6,str7,ID);

else
    fprintf(FD31,'%s\n',str1);
    fprintf(FD31,'%s %s %d %d',str2,file_name,1,length(Seq));
end

fclose(FD31); 
fclose(FD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%getjson%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['Print GI visualization\n' ]);

outjson=strcat(file_name,'_part','.html');
FDjs=fopen(outjson,'w');
fprintf(FDjs,'%s\n','var a={');
fprintf(FDjs,'%s\n','"name":"flare",');
fprintf(FDjs,'%s\n',' "children":[');

for i=1:length(a)-1


if isempty(b{i})
    fprintf(FDjs,'%s\n','{');
    fprintf(FDjs,'%s%s%s%s%s\n',' "name":"',num2str(i),'","size":', num2str(a(i)),',"isco":1,"co":"d3.rgb(209,209,209)"');
    fprintf(FDjs,'%s\n','},');
else
    fprintf(FDjs,'%s\n','{');
    fprintf(FDjs,'%s%s%s%s%s%s%s%s%s\n',' "name":"',num2str(i),'","ss":',num2str(allcanshu(i,1)),',"ee":',num2str(allcanshu(i,2)),',"ll":',num2str(allcanshu(i,4)),',"cg":',num2str(allcanshu(i,5)),',"isco":1,"co":"d3.rgb(228,95,71)",');
    fprintf(FDjs,'%s\n','  "children":[');
    for j=1:length(b{i})-1
    fprintf(FDjs,'%s%s%s%s%s%s%s%s%s%s%s%s%s\n',' {"name":"',strcat(num2str(i),num2str(j)),'","size":',num2str(b{i,1}(j)),',"isco":0,"co":"0","start":',num2str(b{i,4}(2*j-1)),',"end":"',num2str(b{i,4}(2*j)),'","possible":',num2str(b{i,2}(j)),',"ccgg":',num2str(d{i,1}(j)),',"L":',num2str(L),'},');
    end
    fprintf(FDjs,'%s%s%s%s%s%s%s%s%s%s%s%s%s\n',' {"name":"',strcat(num2str(i),num2str(length(b{i}))),'","size":',num2str(b{i,1}(length(b{i}))),',"isco":0,"co":"0","start":',num2str(b{i,4}(2*length(b{i})-1)),',"end":"',num2str(b{i,4}(2*length(b{i}))),'","possible":',num2str(b{i,2}(length(b{i}))),',"ccgg":',num2str(d{i,1}(length(b{i}))),',"L":',num2str(L),'},');
    fprintf(FDjs,'%s%s%s%s%s\n',' {"name":"',strcat(num2str(i),num2str(length(b{i})+1)),'","size":',num2str(int32(sum(b{i,1}(:))*0.02)),',"isco":1,"co":"d3.rgb(255,255,255)"}');
	fprintf(FDjs,'%s\n',']');
     fprintf(FDjs,'%s\n',' },');
end

end
%%%%%last one%%%%%%%%%%%%%%%%%%%%%
i=length(a);

if isempty(b{i})
    fprintf(FDjs,'%s\n','{');
    fprintf(FDjs,'%s%s%s%s%s\n',' "name":"',num2str(i),'","size":', num2str(a(i)),',"isco":1,"co":"d3.rgb(209,209,209)"');
    fprintf(FDjs,'%s\n','},');
	fprintf(FDjs,'%s%s%s%s%s\n',' {"name":"',strcat(num2str(i),num2str(length(b{i})+1)),'","size":',num2str(int32(length(Seq)*0.02)),',"isco":1,"co":"d3.rgb(255,255,255)"}');
else
    fprintf(FDjs,'%s\n','{');
    fprintf(FDjs,'%s%s%s%s%s%s%s%s%s\n',' "name":"',num2str(i),'","ss":',num2str(allcanshu(i,1)),',"ee":',num2str(allcanshu(i,2)),',"ll":',num2str(allcanshu(i,4)),',"cg":',num2str(allcanshu(i,5)),',"isco":1,"co":"d3.rgb(228,95,71)",');
    fprintf(FDjs,'%s\n','  "children":[');
    for j=1:length(b{i})-1
    fprintf(FDjs,'%s%s%s%s%s%s%s%s%s%s%s%s%s\n',' {"name":"',strcat(num2str(i),num2str(j)),'","size":',num2str(b{i,1}(j)),',"isco":0,"co":"0","start":',num2str(b{i,4}(2*j-1)),',"end":"',num2str(b{i,4}(2*j)),'","possible":',num2str(b{i,2}(j)),',"ccgg":',num2str(d{i,1}(j)),',"L":',num2str(L),'},');
    end
    fprintf(FDjs,'%s%s%s%s%s%s%s%s%s%s%s%s%s\n',' {"name":"',strcat(num2str(i),num2str(length(b{i}))),'","size":',num2str(b{i,1}(length(b{i}))),',"isco":0,"co":"0","start":',num2str(b{i,4}(2*length(b{i})-1)),',"end":"',num2str(b{i,4}(2*length(b{i}))),'","possible":',num2str(b{i,2}(length(b{i}))),',"ccgg":',num2str(d{i,1}(length(b{i}))),',"L":',num2str(L),'},');
    fprintf(FDjs,'%s%s%s%s%s\n',' {"name":"',strcat(num2str(i),num2str(length(b{i})+1)),'","size":',num2str(int32(sum(b{i,1}(:))*0.02)),',"isco":1,"co":"d3.rgb(255,255,255)"}');
	
	fprintf(FDjs,'%s\n',']');
     fprintf(FDjs,'%s\n',' }');
end
fprintf(FDjs,'%s\n',']');
fprintf(FDjs,'%s\n','};');

fclose(FDjs);
 
FD4=fopen(out4,'w');

fid1=fopen('part1.html','r');
fid2=fopen('part2.html','r');
FDjs=fopen(outjson,'r');
Data1=fread(fid1);
Data2=fread(fid2);
Datajs=fread(FDjs);
fwrite(FD4,Data1);
fwrite(FD4,Datajs);
fwrite(FD4,Data2);
fclose(fid1);
fclose(fid2);
fclose(FDjs);
fclose (FD4);
delete('*_part.html')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  fclose(FD2); 






