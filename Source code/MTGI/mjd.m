function score_js=mjd(Seq,position)

NU='ATCG';
seq1=Seq(1:position);
seq2=Seq(position+1:length(Seq));
fren=zeros(3,4^3);
tran=zeros(3,4^3);

% the 2-order Markov model of Seq
    fren_s=zeros(4,4,4); 
    tran_s=zeros(4,4,4);
    for k1=1:4 % transition probablity
        for k2=1:4
               for k3=1:4
                   a=[];b=[];
                   a=[NU(k1) NU(k2) NU(k3)];
                   b=findstr(Seq,a);
                   tran_s(k1,k2,k3)=length(b)+1;
               end
         end
    end
    fren_s=tran_s/(length(Seq)+64-2);
    fren(1,:)=reshape(fren_s,1,64);
    for k1=1:4
           for k2=1:4
               if sum(tran_s(k1,k2,:))~=0
                  tran_s(k1,k2,:)=tran_s(k1,k2,:)/sum(tran_s(k1,k2,:));
               end
           end       
    end
    tran(1,:)=reshape(tran_s,1,64);
    
    % the 2-order Markov model of Seq1
    fren_s1=zeros(4,4,4); 
    tran_s1=zeros(4,4,4);
    for k1=1:4 % transition probablity
        for k2=1:4
               for k3=1:4
                   a=[];b=[];
                   a=[NU(k1) NU(k2) NU(k3)];
                   b=findstr(seq1,a);
                   tran_s1(k1,k2,k3)=length(b)+1;
               end
         end
    end
    fren_s1=tran_s1/(length(seq1)+64-2);
    fren(2,:)=reshape(fren_s1,1,64);
    for k1=1:4
           for k2=1:4
               if sum(tran_s1(k1,k2,:))~=0
                  tran_s1(k1,k2,:)=tran_s1(k1,k2,:)/sum(tran_s1(k1,k2,:));
               end
           end       
    end
   tran(2,:)=reshape(tran_s1,1,64);
    
% the 2-order Markov model of Seq2
    fren_s2=zeros(4,4,4); 
    tran_s2=zeros(4,4,4);
    for k1=1:4 % transition probablity
        for k2=1:4
               for k3=1:4
                   a=[];b=[];
                   a=[NU(k1) NU(k2) NU(k3)];
                   b=findstr(seq2,a);
                   tran_s2(k1,k2,k3)=length(b)+1;
               end
         end
    end
    fren_s2=tran_s2/(length(seq2)+64-2);
    fren(3,:)=reshape(fren_s2,1,64);
    for k1=1:4
           for k2=1:4
               if sum(tran_s2(k1,k2,:))~=0
                  tran_s2(k1,k2,:)=tran_s2(k1,k2,:)/sum(tran_s2(k1,k2,:));
               end
           end       
    end
    tran(3,:)=reshape(tran_s2,1,64);
 
 % calculate JS divergence of two segmentation of Seq
 score_js=-sum(fren(1,:).*log(tran(1,:)))+position/length(Seq)*sum(fren(2,:).*log(tran(2,:)))+(length(Seq)-position)/length(Seq)*sum(fren(3,:).*log(tran(3,:)));
end