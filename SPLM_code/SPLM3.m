%%
clear;clc;
na=load('ATPB647+TUB660+TOM5.5_SR10nm_5ms 2018 April 22 14_40_53_TR');
bs=load('ATPB647+TUB660+TOM5.5_SR10nm_5ms 2018 April 22 14_40_53_spe');
na0=load('ATPB647+TUB660+TOM5.5_SR10nm_5ms 2018 April 22 14_40_53_TR0');
img=(imread('ATPB647+TUB660+TOM5.5_SR10nm_5ms 2018 April 22 14_40_53.tif'));
%%

bs(isnan(bs))=0;
s=sum(bs,2);
s=s/100*1/0.925;

%Flitering based on spectral intensity
a1=s>=0;
a2=s<=1e5;
a=logical(a1.*a2);
s=s(a);
spe=bs(a,:);
no=na0(a,3:4); % coordinate

wl=[650:2:775];
c=[];
for n=1:size(spe,1)
    c(n)=sum(wl.*(spe(n,:).^1))./sum(spe(n,:).^1);
end

figure(1);
clf
fy = [650:25:775];%for xtick in plotting

subplot 231
imshow(img);

subplot 232
scatter(c,s,4,'bo','filled')
axis([650 750 0 1e4])
ylabel('Intenisty')
xlabel('Wavelength (nm)')
set(gca,'TickLength',[0.02,0.01]);
box on

subplot 233
[nb,xb]=hist(c,[600:1:800]); 
b=bar(xb,nb);
set(b,'facecolor',[1,1,1]*0.5);

set(gca,'TickLength',[0.02,0.01]);
xlabel('Wavelength (nm)')
ylabel('Count')
xlim([650 750])
std(c)

%%
%SELECT ROIS IN ORDER
new_ROIs = 1;

%%Cycle through ROIs and obtain mask for each if needed
if new_ROIs == 1
     z=figure(10);
%     h_im = imagesc(I{13}); %Show strongest signal image
    h_im = imshow(img);
    colormap gray
    text(0,-30,'Left click, adjust Rectangular ROI and left click twice, and wait for result.',...
        'FontSize',[10],'Color', 'b');
    
    roi = imrect;  % get a region R inside a rectangle
    wait(roi);
    ImgMaskroi = createMask(roi,h_im);
    save('ImgMaskroi.mat','ImgMaskroi');%save to pwd
   close(z);
    roic=ans*10;
else
  % load saved ROI masks if not making new ones
    disp('Reading from pwd; if this does not work set new_ROIs = 1 above')
    load('ImgMaskroi.mat')
end
    
%%
out_new=double(img(:,:,1)).*ImgMaskroi;
figure(1)
subplot 234
imshow(out_new*100);  
axis image;

row_start =  roic(1);
row_end = roic(1)+roic(3);
col_start =  roic(2);
col_end = roic(2)+roic(4);

a1=no(:,1)>col_start;
a2=no(:,1)<col_end;
a3=no(:,2)>row_start;
a4=no(:,2)<row_end;
a=logical(a1.*a2.*a3.*a4);
no_new=no(a,:);
s_new=s(a);
c_new=c(a);
fy = [650:25:775];%for xtick in plotting


subplot 235
scatter(c_new,s_new,4,'bo','filled')
axis([650 750 0 2e3])
ylabel('Intenisty')
xlabel('Wavelength (nm)')
set(gca,'TickLength',[0.02,0.01]);
box on

subplot 236
[nb,xb]=hist(c_new,[600:1:800]); 
b=bar(xb,nb);
set(b,'facecolor',[1,1,1]*0.5);

set(gca,'TickLength',[0.02,0.01]);
xlabel('Wavelength (nm)')
ylabel('Count')
xlim([650 750])
std(c_new)

