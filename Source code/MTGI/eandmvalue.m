function ScoringList=eandmvalue(data,main_data,eye,ChosenTetra)

ScoringList=zeros(size(data,1),1);
if nargin<4
    ChosenTetra=1:size(data,2);
end
m=[];sig=[];
m=mean(main_data);
sig=std(main_data);
Num_feature=length(ChosenTetra);
for l=1:Num_feature
    slide_value=zeros(2*eye+1,2*eye+size(data,1));
    for zz=1:2*eye+1
        slide_value(zz,:)=[zeros(1,2*eye+1-zz) data(:,ChosenTetra(l))' zeros(1,zz-1)];
    end
    average_eye=[];variance_eye=[];
    average_eye=mean(slide_value)';variance_eye=var(slide_value)';
    average_eye(1:eye)=[];
    average_eye(size(data,1)+1:eye+size(data,1))=[];
    sp=[];
    variance_eye(1:eye)=[];
    variance_eye(size(data,1)+1:eye+size(data,1))=[];
    sp=sqrt((2*eye.*variance_eye+(size(main_data,1)-1)*sig(ChosenTetra(l))^2)/(2*eye+size(main_data,1)-1));
    tem_vec=[];
    tem_vec=(average_eye-m(ChosenTetra(l)))./(sp*sqrt((2*eye+1+size(main_data,1))/(size(main_data,1)*(2*eye+1))));
    score_tem=-log(1-tcdf(abs(tem_vec),size(main_data,1)+2*eye-1));
    % Impute inf data using nearest-neighbor value
    Score_impute=[];
    Score_impute=score_tem;
    wid_value=10;
    position_inf=find(score_tem==inf);
    Score_impute(position_inf)=0;
    for num_position_inf=1:length(position_inf);
        start_window=max(position_inf(num_position_inf)-wid_value,1);
        end_window=min(position_inf(num_position_inf)+wid_value,length(score_tem));
        Score_impute(position_inf(num_position_inf))=mean(Score_impute(start_window:end_window));
    end
    ScoringList=ScoringList+Score_impute;
end
end