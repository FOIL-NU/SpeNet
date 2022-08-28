clear
clc
clf

N=100000;

sect=7;
r0=1;
d=0.1;
shift=0.08;
center=3+4*1i;

NA=3*rand(N,1)+6*1i*rand(N,1);

NA(real(NA)<0)=[];
NA(real(NA)>3)=[];
NA(imag(NA)<0)=[];
NA(imag(NA)>6)=[];

NA(abs(NA)<0.2)=0;

NA2=[];
for n=1:sect
    l1=abs(NA-center+(2*n-2)*shift*1i)>=r0+d*(n*2-2);
    l2=abs(NA-center+(2*n-1)*shift*1i)<=r0+d*(n*2-1);
    l=logical(l1.*l2);
    NA2=[NA2;[real(NA(l)),imag(NA(l)),ones(sum(l),1)*n]];
end

scatter(NA2(:,1),NA2(:,2),1,sect+1-NA2(:,3),'filled')
colormap(jet)
axis equal
axis([0 3 0 6])
%%
NA2=NA2(randperm(size(NA2,1)),:);
save('NA_c.txt','NA2','-ascii')
%%
%draw
clf
clc

NA=load('NA_c.txt');
d=0.002;
shift=0;
color=jet(7);

out=zeros((3+shift)/d+1,(6+shift)/d+1,3);

for k=1:3
    for n=1:7
        NS=NA(NA(:,3)==n,:);
        out(:,:,k)=out(:,:,k)...
            +color(n,k)*hist3([NS(:,1)+shift,NS(:,2)+shift],{0:d:3+shift 0:d:6+shift});
    end
end

imshow(out)

%%
out2=zeros((3+shift)/d+1,(6+shift)/d+1,3);
gs=11;
g=fspecial('gaussian', gs, 2);

for k=1:3
    out2(:,:,k)=conv2(out(:,:,k),g,'same');
end

% out2=out2([1:600],[1:1200],:);
imshow(out2/max(out2(:))*3)
%%
imwrite(out2/max(out2(:))*2,'Img_2.tif')
%%
addpath ('E:\Matlab_code\SPLM_code'); 
clear

na=load('na.dat');
bs=load('spe.dat');
%%

bs(isnan(bs))=0;
s=sum(bs,2);
s=s/100*3/0.95;

%Flitering based on spectral intensity
a1=s>=400;
a2=s<=10000;
% a3=sum(bs(:,1:10),2)<=30;
% a4=sum(bs(:,end-30:end),2)<=500;
% a=logical(a1.*a2.*a3.*a4);
a=logical(a1.*a2);
s=s(a);
spe=bs(a,:);
no=na(a,1:2); % coordinate

wl=[525:1:700];
c=[];
for n=1:size(spe,1)
    c(n)=sum(wl.*(spe(n,:).^1))/sum(spe(n,:).^1);
end

%%
a1=c>710-30;
a2=c<710+30;
a=logical(a1.*a2);
c2=c(a);
no2=no(a,:);
c2=c2-660;
c2=round(((c2+20.5)));

cr = [];
cmap = jet (101);
% cmap = cmap ([10:20],:);
for n = 1:size (c2,2)
    rgb = cmap (c2(n),:);
    cr = [cr;rgb];
end

scatter (no2(:,1),no2(:,2),1,cr);
%%
temp = zeros (30,30,3);
temp (2,2) = 1;
m = 100;
temp (:,:,1) = cmap (m,1);
temp (:,:,2) = cmap (m,2);
temp (:,:,3) = cmap (m,3);
imshow (temp)

%%

x = round (no2 (:,1)*10);
y = round (no2 (:,2)*10);
%%
out =  zeros (3500,3500,3);
out2 = out;
% x = round (x/100);
% y = round (y/100);
for n = 1:length (no2)
    for k=1:3
    out (x(n),y(n),k)=cr (n,k);
    end
end

gs=50;
g=fspecial('gaussian', gs, 5);

% out2=zeros (2e3,2e2,3);

for k=1:3
    out2(:,:,k)=conv2(out(:,:,k),g,'same');
end
rgb1=uint8 (out2);
%%
image = imshow (out2*100);
%%
imwrite(out2*2,'Spatial_Tub+TOM20_10minbleach022_mito2.tif')

%%
NA=[no2,c2'];
color=jet(100);


d=0.002;
cbin = 100;
  
% color=color([1:90]+10,:);

out=zeros(40000+1,20000+1,3);

for k=1:3
    for n=1:cbin
        NS=NA(NA(:,3)==n,:);
        out(:,:,k)=out(:,:,k)...
            +color(n,k)*hist3([NS(:,1),NS(:,2)],{0:10:400000,0:10:200000});
    end
end
%%
imshow(out)
%%
out2=zeros(4e3+1,2e3+1,3);

gs=8;
g=fspecial('gaussian', gs, 5);

for k=1:3
    out2(:,:,k)=conv2(out(:,:,k),g,'same');
end

% out2=out2([1:600],[1:1200],:);
imshow(out2/max(out2(:))*3)
%%
imwrite(out2/max(out2(:))*2,'clathrino-3.tif')

%%
imagesc([1:1:1000])
colormap(jet(100))

