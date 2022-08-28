%% normalize data to [0,1]
a=load('680 3_R0.95.mat');b=load('680 2_R0.95.mat');
an680=cat(1,a.spe,b.spe);
%%
clear
an647=load ('647_140mW_20ms001_specimg_full.mat');an660=load('660_140mW_20ms003_specimg_full.mat');an680=load('680_140mW_20ms_006_specimg_full.mat');
%% intensity filter
s680=an680.s;
an680=s680.raw_speimage;
i680=s680.i;
an680=an680(i680>2000,:,:);
m = size (an680,1);
an680_n = zeros (m,7,201);
for n = 1: m
    temp = squeeze (an680(n,:,:));
    temp1 = temp./max(temp(:));
    an680_n (n,:,:) = temp1;
end
s647=an647.s;
an647=s647.raw_speimage;
i647=s647.i;
an647=an647(i647>2000,:,:);
m = size (an647,1);
an647_n = zeros (m,7,201);
for n = 1: m
    temp = squeeze (an647(n,:,:));
    temp1 = temp./max(temp(:));
    an647_n (n,:,:) = temp1;
end
s660=an660.s;
an660=s660.raw_speimage;
i660=s660.i;
an660=an660(i660>2000,:,:);

m = size (an660,1);
an660_n = zeros (m,7,201);
for n = 1: m
    temp = squeeze (an660(n,:,:));
    temp1 = temp./max(temp(:));
    an660_n (n,:,:) = temp1;
end
%% coefficient filter
R647=an647.s.R647;
sub_speimage=an647.s.sub_speimage;
i=an647.s.i;
m = length (R647);
sub_speimage_n = zeros (m,7,201);
for n = 1: m
    temp = squeeze (sub_speimage(n,:,:));
    temp1 = temp./max(temp(:));
    sub_speimage_n (n,:,:) = temp1;
end
an647_n=sub_speimage_n(R647>0.7,:,:);
% an_noise1_n=raw_speimage_n(R647<=0.7,:,:);
%%
R660=an660.s.R660;
raw_speimage=an660.s.raw_speimage;
i=an660.s.i;
m = length (R660);
raw_speimage_n = zeros (m,7,201);
for n = 1: m
    temp = squeeze (raw_speimage(n,:,:));
    temp1 = temp./max(temp(:));
    raw_speimage_n (n,:,:) = temp1;
end
an660_n=raw_speimage_n(R660>0.8,:,:);
an_noise2_n=raw_speimage_n(R660<=0.8,:,:);
%%
R680=an680.s.R680;
raw_speimage=an680.s.raw_speimage;
i=an680.s.i;
m = length (R680);
raw_speimage_n = zeros (m,7,201);
for n = 1: m
    temp = squeeze (raw_speimage(n,:,:));
    temp1 = temp./max(temp(:));
    raw_speimage_n (n,:,:) = temp1;
end
an680_n=raw_speimage_n(R680>0.8,:,:);
an_noise3_n=raw_speimage_n(R680<=0.8,:,:);
%%
an647r= permute (an647_n,[2,3,1]);
an660r= permute (an660_n,[2,3,1]);
an680r= permute (an680_n,[2,3,1]);
% an647r=an680r; % for 660-680 validation only, 647 are now 680 data
clearvars -except an647r an660r an680r;
% clearvars -except an647r an660r;
XTrain = cat (3,an647r,an660r,an680r);
% XTrain = cat (3,an647r,an660r);
% XTrain = cat (3,an660r,an680r);
XTrain = permute (cat(4,XTrain,XTrain),[1,2,4,3]);
XTrain = XTrain (:,:,1,:);
a1 = zeros (size(an647r,3),1);
a2 = ones (size(an660r,3),1);
a3 = (ones (size(an680r,3),1))*2;
YTrain = cat (1,a1,a2,a3);
% YTrain = cat (1,a2,a3);
%%% shuffle
IDX=randperm(length(XTrain));
XTrain=XTrain(:,:,:,IDX);
YTrain= categorical(YTrain(IDX));

%%%
idx = randperm(size(XTrain,4),round(size(XTrain,4)*0.2)); % 20% for validation
XValidation = XTrain(:,:,:,idx);
XTrain(:,:,:,idx) = [];
YValidation = YTrain(idx);
YTrain(idx) = [];

%%%
layers = [
    imageInputLayer([7 201 1])
    
    convolution2dLayer(5,8,'Padding','same')
    batchNormalizationLayer
    reluLayer
%     dropoutLayer
    
%     maxPooling2dLayer(2,'Stride',2)
%     
%     convolution2dLayer(3,16,'Padding','same')
%     batchNormalizationLayer
%     reluLayer   
    
    maxPooling2dLayer(3,'Stride',2)
    
    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer   
%     dropoutLayer

    fullyConnectedLayer(3)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm', ...
    'MaxEpochs',5, ...
    'ValidationData',{XValidation,YValidation}, ...
    'ValidationFrequency',30, ...
    'ExecutionEnvironment','gpu',...
    'Verbose',false, ...
    'Shuffle','once',...
    'Plots','training-progress');
%%%% train network
net = trainNetwork(XTrain,YTrain,layers,options);
%% save net
save ('CNN_2021_07_12.mat','net');

%%
clf
[YPred, scores] = classify (net,XValidation);
accuracy=sum(YPred==YValidation)/length(YValidation);
C = confusionmat(YValidation,YPred);
confusionchart(C)
%% systematic validation
a = load ('680_140mW_20ms_006_specimg_full.mat');
b = load ('680_140mW_20ms_006_specimg_full.mat');
photon647=[a.s.i;b.s.i];%;c.s.i;d.s.i];
centroid647=[a.s.c;b.s.c];%;c.s.c;d.s.c];
R647=[a.s.R647;b.s.R647];%;c.s.R647;d.s.R647];
R660=[a.s.R660;b.s.R660];%c.s.R660;d.s.R660];
R680=[a.s.R680;b.s.R680];%c.s.R680;d.s.R680];
spe647=cat(1,a.s.raw_speimage,b.s.raw_speimage);%c.s.raw_speimage,d.s.raw_speimage);
%%
acc647=[];
acc660=[];
acc680=[];
uti=[];
% for n = [-1,0.5:0.05:0.95]
 spe647_f=spe647(R647>-0.8,:,:);
a = permute (cat(4,spe647_f,spe647_f),[2,3,4,1]);
a = a (:,:,1,:);
[YPred, scores]=classify (net,a);
v=categorical(zeros (length(spe647_f),1));
v1=categorical(ones(length(spe647_f),1));
v2=categorical(2.*ones(length(spe647_f),1));
acc647=[acc647;sum(YPred==v)/length(spe647_f)];
acc660=[acc660;sum(YPred==v1)/length(spe647_f)];
acc680=[acc680;sum(YPred==v2)/length(spe647_f)];
uti=[uti;length(v)./length(spe647)];
% end
acc=cat(2,acc647,acc660,acc680,uti);