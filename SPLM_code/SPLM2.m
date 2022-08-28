%%
clear
fname = 'Tubuline_Alexa647_SR3nm_2018 April 26 19_03_59';
save1 = sprintf('%s_TR',fname);
save2 = sprintf('%s_TR0',fname);
save3 = sprintf('%s_spe',fname);
save4 = sprintf('%s_bb3',fname);

na=load(save1);
bs=load(save3);
na0=load(save2);
bb3=load (save4)';
%%
bs(isnan(bs))=0;
s=sum(bb3,2)/100*4.6/0.95;
% s=na0(:,6);

%Flitering based on spectral intensity
a1=s>=1500;
a2=s<=2e6;
a=logical(a1.*a2);
s=s(a);
spe=bs(a,:);
no3=na0(a,:); % coordinate
no=na0(a,3:4); % coordinate

wl=[600:1:800];
c=[];
for n=1:size(spe,1)
    c(n)=sum(wl.*(spe(n,:).^1))./sum(spe(n,:).^1);
end


%%
tic;
options = fitoptions('gauss1', 'Lower', [0 -inf 0 ],'Upper',[inf inf 30]);
for i=1:length(spe)
    warning off;
    f = fit(wl',spe(i,:)','gauss1',options);
    cx2(i,:) = f.b1;
end
toc;
% 
% tic;
% options = optimset('Display','off','TolFun',4e-16,'LargeScale','off');
% for i=1:length(spe)
% input = spe(i,:);
% data=input.*(input>0);
% x1 = [1:length(data)];
% cx1 = sum(wl.*(spe(i,:).^1))/sum(spe(i,:).^1);
% sx1 = sqrt(sum((abs((wl-cx1).^2).*data)/sum(data)));
% peak1 = max(max(data));
% 
% para = [cx1,sx1,peak1];
% fp1D = fminunc(@fitgaussian1D,para,options,data,x1);
% cx2(i) = fp1D(1)';
% sx2 = fp1D(2);
% peak2 = fp1D(3);
% 
% end
% toc;
c=cx2';
% save('c_fitting.txt','c','-ascii')
% c=load('c_fitting.txt');
%%
for nn = 705
wl1=nn;
wl2=nn+6;
wln=wl2-wl1;
a1=c>wl1;
a2=c<wl2;
a=logical(a1.*a2);
c2=c(a);
no2=no(a,:);
no4=no3(a,:);
spe1=spe(a,:);
c2=round(c2);
c2=c2-min(c2);
NA=[no2,c2'];

M = 1.05; % magnification factor of window
xwindow=round(max(na0(:,3))/10*M);
ywindow=round(max(na0(:,4))/10*M);
out=zeros(xwindow+1,ywindow+1,3);

color=jet(wln*20);
color=color([1:max(c2)]+wln*20-max(c2),:);

for k=1:3
    for n=1:max(c2)
        NS=NA(NA(:,3)==n,:);
        out(:,:,k)=out(:,:,k)+color(n,k)*hist3([NS(:,1),NS(:,2)],{0:10:xwindow*10,0 :10:ywindow*10});
    end
end

out2=zeros(xwindow+1,ywindow+1,3);
g=fspecial('gaussian', 7, 3);

for k=1:3
    out2(:,:,k)=conv2(out(:,:,k),g,'same');
end

figure(nn)
% title (n);
imshow(out2/max(out2(:))*20);
drawnow;
pause (1);
end
%%
wl1=;
wl2=712;
wln=wl2-wl1;
a1=c>wl1;
a2=c<wl2;
a=logical(a1.*a2);
c2=c(a);
no2=no(a,:);
c2=round(c2);
c2=c2-min(c2);
NA=[no2,c2'];


out3=zeros(xwindow+1,ywindow+1,3);
color=jet(wln*20);
color=color([1:max(c2)]+wln*20-8*max(c2),:);

for k=1:3
    for n=1:max(c2)
        NS=NA(NA(:,3)==n,:);
        out3(:,:,k)=out3(:,:,k)+color(n,k)*hist3([NS(:,1),NS(:,2)],{0:10:xwindow*10,0 :10:ywindow*10});
    end
end

g=fspecial('gaussian', 7, 3);
out4=zeros(xwindow+1,ywindow+1,3);

for k=1:3
    out4(:,:,k)=conv2(out3(:,:,k),g,'same');
end

imshow(out4/max(out4(:))*40);
%%

out5 = out2/max(out2(:))*20;
out6 = out4/max(out4(:))*30;
C = imfuse(out5,out6,'falsecolor','Scaling','joint','ColorChannels',[2 1 2]);
% C = imfuse(out5,out6,'blend');
figure(2); 
imshow(out6);hold on;
%%
fnamein=sprintf('683-688-400.tif',fname);
imwrite(out2,fnamein);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Ver. 2 %%%%%%%%%%%%%%%%%%%%%%%

bs(isnan(bs))=0;
s=sum(bs,2);
s=s/100*3/0.925;

%Flitering based on spectral intensity
a1=s>=200;
a2=s<=10000;
a=logical(a1.*a2);
s=s(a);
spe=bs(a,:);
no=na0(a,3:4); % coordinate

wl=[650:2:775];
c=[];
for n=1:size(spe,1)
    c(n)=sum(wl.*(spe(n,:).^1))/sum(spe(n,:).^1);
end

%%
options = fitoptions('gauss1', 'Lower', [0 -inf 0 ],'Upper',[inf inf 30]);
for i=1:length(spe)
    warning off;
    f = fit(wl',spe(i,:)','gauss1',options);
    cx2(i,:) = f.b1;
end
c=cx2';
%%
ref1=690;
ref2=30;
a1=c>ref1-ref2;
a2=c<ref1+ref2;
a=logical(a1.*a2);
c2=c(a);
no2=no(a,:);
c2=c2-(ref1-ref2);
c2=round(((c2/max(c2)))*160+0.5);

cr = [];
cmap = jet (max(c2));
% cmap = cmap ([10:20],:);
for n = 1:size (c2,2)
    rgb = cmap (c2(n),:);
    cr = [cr;rgb];
end

scatter (no2(:,1),no2(:,2),1,cr);
%%
x = round (no2 (:,1)/160*10+0.5);
y = round (no2 (:,2)/160*10+0.5);
out =  zeros (max(x),max(y),3);

out2 = out;
for n = 1:length (no2)
    for k=1:3
    out (x(n),y(n),k)=cr (n,k);
    end
end

g=fspecial('gaussian', 5, 5);

for k=1:3
    out2(:,:,k)=conv2(out(:,:,k),g,'same');
end
figure(4);
image = imshow (out2*5);