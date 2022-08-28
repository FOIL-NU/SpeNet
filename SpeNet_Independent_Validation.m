%% input data by clicking and filter validation data by coefficient

R1=0;R2=0;R3=0;
s6471=s647_v(R647_v(:,1)>R1 | R647_v(:,2)>R2 | R647_v(:,3)>R3,:,:);
m = length (s6471);
s647_n = zeros (m,7,201);
for n = 1: m
    temp = squeeze (s6471(n,:,:));
    temp1 = temp./max(temp(:));
    s647_n (n,:,:) = temp1;
end
s6601=s660_v(R660_v(:,1)>R1 | R660_v(:,2)>R2 | R660_v(:,3)>R3,:,:);
m = length (s6601);

s660_n = zeros (m,7,201);
for n = 1: m
    temp = squeeze (s6601(n,:,:));
    temp1 = temp./max(temp(:));
    s660_n (n,:,:) = temp1;
end

s6801=s680_v(R680_v(:,1)>R1 | R680_v(:,2)>R2 | R680_v(:,3)>R3,:,:);
m = length (s6801);
s680_n = zeros (m,7,201);
for n = 1: m
    temp = squeeze (s6801(n,:,:));
    temp1 = temp./max(temp(:));
    s680_n (n,:,:) = temp1;
end  
%%% combine data for network classifying
a1 = categorical(zeros (length(s647_n),1));
a2 = categorical(ones (length(s660_n),1));
a3 = categorical((ones (length(s680_n),1))*2);
raw_speimage_n=cat(1,s647_n,s660_n,s680_n);
spe=cat(4,raw_speimage_n,raw_speimage_n);
spe=permute(spe(:,:,:,1),[2,3,4,1]);
YValidation=cat(1,a1,a2,a3);
%%% simple classifying
clf
[YPred,scores]=classify(net,spe);
C = confusionmat(YValidation,YPred); confusionchart(C,'Normalization','column-normalized')
% a647_660=C(1,2)/C(1,1);a660_647=C(2,1)/C(2,2);a660_680=C(2,3)/C(2,2);a680_660=C(3,2)/C(3,3);
% accuracy=sum(YPred==YValidation)/length(YValidation);
%%% pick classifying results with prediction score filter
C_w=[];
for f=0.35:0.05:0.95
[m,i]=max(scores,[],2);
idx=m>f;i1=categorical(i(idx)-1);
YPred_f=YPred(idx);
YValidation_f=YValidation(idx);
accuracy=sum(YPred_f==YValidation_f)/length(YValidation_f);
C = confusionmat(YValidation_f,YPred_f); 
% C (:,1)=C(:,1)/sum(YValidation_f==categorical(0))*100;C(:,2)=C(:,2)/sum(YValidation_f==categorical(1))*100;C(:,3)=C(:,3)/sum(YValidation_f==categorical(2))*100;
confusionchart(C,'Normalization','row-normalized');
C_w=[C_w;C];
end
%%%
A647=C_w(1:3:end,:);
A660=C_w(2:3:end,:);
A680=C_w(3:3:end,:);
AU=[A647,A660,A680];