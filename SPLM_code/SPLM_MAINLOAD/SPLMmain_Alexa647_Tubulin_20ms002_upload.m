%Calibration data
yi2  = [139, 154];
xi=[531,617.8];
res = (peak(end)-peak(1))/(yi2(end)-yi2(1)) % res = 5.7867

%%
% Load spectrum file
clear
clc;
addpath('E:\Matlab_code\SPLM_code\SPLM_MAINLOAD');
fnamein=sprintf('3color_10ms004.nd2');
A = SPLMload(fnamein,'nd2','double');
A = rot90(A,3);
figure(1);imagesc(squeeze(mean(A,3))); axis off; axis image; axis image; colormap(gray)
A=double(squeeze(A));
A= A - 200;
%%

clf
for n=1:500
    imagesc(A(:,:,n));
    title(n)
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    colormap(gray)
    set(gca,'color','red','Position',[0.05 0.1 0.9 0.8]);
    drawnow
%     pause(0.01)
end

x1=10;
x2=90;
%%
%Cropping image for ThunderStorm reconstruction
fnameout=sprintf('3color_10ms004.dat');
img=A(:,x1:x2,:);
fid1 = fopen(fnameout, 'w') ;
fwrite(fid1, img, 'uint16') ;
fclose(fid1);


imagesc(img(:,:,1));axis equal;
%% Optional
%Calculate gobal background
% 
back=[];
for n=1:size(A,1)
    w=squeeze(A(n,x1:x2,:));
    w=mean(w,1);
    a=w<500;
    back(n,:)=mean(A(n,:,a),3);
end
colormap gray
imagesc(back)
%%

%Load reconstruction result from Thunderstorm
TR0=importdata('3color_10ms004.csv');
TR0=TR0.data;

TR=TR0(:,[3,4,2,6]);
PixelSize=160;
TR(:,1)=round(TR(:,1)/PixelSize+0.5);
TR(:,2)=round(TR(:,2)/PixelSize+0.5);

%Cropping for further process
w=4;
ax_min=TR(:,1)>w+1;
ax_max=TR(:,1)<max(TR(:,1))-w-1;
ay_min=TR(:,2)>w+1;
ay_max=TR(:,2)<max(TR(:,2))-w-1;
az_min=TR(:,3)>=w+1;
az_max=TR(:,3)<=max(TR(:,3))-w;
ai=TR(:,4)>=0;
a=ax_min.*ax_max.*ay_min.*ay_max.*az_min.*az_max.*ai;
a=logical(a);
TR=TR(a,:);
TR0=TR0(a,:);

out=hist3([TR0(:,3),TR0(:,4)],{0:10:max(TR0(:,3))*1.05,0:10:max(TR0(:,4)*1.05)});
figure (2);
imshow(out)

%%

%Spectral calibration data (three wavelengths)
yi=[102,118];
xi=[570,671];
fy = [625:25:775];%for xtick in plotting
fy2= [625:2:775];%for interpolation
coeffs=polyfit(xi, yi, 1);
fx = polyval(coeffs, fy);
fx2 = polyval(coeffs, fy2);


x0=0;%slit position
xc=x1-1;%cropping x-position

y0=0;%Y shift

% input
% TR=TRnew;
s=TR(:,[2,1,3,4]);%this order relies on the output format of ThunderStorm
% s=flipud(sortrows(s,4)); %sort by intensity
% bg=permute(bg,[2 1 3]);
s(:,1)=s(:,1)+xc;
s(:,2)=s(:,2);

aa=[];%Blinking image
bb=[];%Spectrum image

as=[];%blinking profile
bs=[];%spectrum profile

bw=1;%blinking half width (actual width of the PSF=2*bw+1)

w=3;%half width of the calculated region
% lbs=[-4:-3,3:4];%for local background subtraction

live=0;%live update: show=1; no show=0 

for n=1: size(s,1)

    % % local background subtraction
    % sz=lbs+s(n,3);
    % a=A([-w:w]+s(n,2),[-w:w]+s(n,1),s(n,3));
    % aa=a-mean(A([-w:w]+s(n,2),[-w:w]+s(n,1),sz),3);
    % b=A([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1),s(n,3));
    % bb=b-mean(A([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1),sz),3);

    % %global background subtraction
    aa=A([-w:w]+s(n,2),[-w:w]+s(n,1),s(n,3))-back([-w:w]+s(n,2),[-w:w]+s(n,1));
    bb=A([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1),s(n,3))-back([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1));

    as(n,:)=mean(aa(w+1-bw:w+1+bw,:));
    bs(n,:)=interp1([round(fx(1)):round(fx(end))],mean(bb(w+1-bw:w+1+bw,:)),fx2);

    % show progress
    if mod(n*10000,size(s,1))<=50
        fprintf('Progress: %d \n',round(n/size(s,1)*100))
    end

    %Live display
    if live==1
        subplot('Position',[0.06 0.55 0.17 0.4])
            imagesc(aa,[0 3000]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);

        subplot('Position',[0.3 0.55 0.65 0.4])
            imagesc(bb,[0 4000])
            set(gca,'XTick',[1:(size(bb,2)-1)/(size(fy,2)-1):size(bb,2)]);
            set(gca,'XTicklabel',fy);
            set(gca,'YTick',[]);
            colormap(gray)

        subplot('Position',[0.06 0.1 0.17 0.35])
            plot(as(n,:),'k');
            % hold on
            % plot(smooth(bs(n,:),20),'r');hold off
            axis([1 size(as,2) -100 3000])
            set(gca,'XTick',[1:(size(bs,2)-1)/(size(fy,2)-1):size(bs,2)]);
            set(gca,'XTicklabel',fy);
        %     set(gca,'YTick',[]); 

        subplot('Position',[0.3 0.1 0.65 0.35])
            plot(bs(n,:),'k');
            % hold on
            % plot(smooth(bs(n,:),20),'r');hold off
            axis([1 size(bs,2) -100 4000])
            set(gca,'XTick',[1:(size(bs,2)-1)/(size(fy,2)-1):size(bs,2)]);
            set(gca,'XTicklabel',fy);
        %     set(gca,'YTick',[]); 
        drawnow
    end
end


%%

%statistic analysis

bs(isnan(bs))=0;
s=sum(bs,2);
s=s/100*4.63/0.95;

%Flitering based on spectral intensity
% a1=s>=500;
% a2=s<=1e5;
% a3=sum(bs(:,1:60),2)<=500;
% a4=sum(bs(:,end-10:end),2)<=1000;
% a=logical(a1.*a2.*a3.*a4);
% % a=logical (a1.*a2);
% s=s(a);
% spe=bs(a,:);

a1=s>=1000;
a2=s<=inf;
% a3=sum(bs(:,1:60),2)<=500;
% a4=sum(bs(:,end-10:end),2)<=1000;
% a=logical(a1.*a2.*a3.*a4);
a=logical (a1.*a2);
s=s(a);
spe=bs(a,:);
% ns=na2(a,:);


%dispersion correction

%replace corr with dispersion calibration data
% corr=meshgrid(linspace(0.7,0.55,451),linspace(1,1,size(spe,1)));
% spe=spe./corr;

 
%Calculate spectral centroid
wl=[625:2:775];
c=[];
for n=1:size(spe,1)
    c(n)=sum(wl.*(spe(n,:).^1))/sum(spe(n,:).^1);
end 
figure(1);
subplot('position',[0.12 0.35 0.28 0.6])

barh(flipud(sum(spe,2)),'EdgeColor','None')
axis([0 max(sum(spe,2)*2) 0.5 size(spe,1)+0.5])
ylabel('Blinking')
xlabel('Intenisty')

subplot('position',[0.46 0.35 0.5 0.6])
imagesc(spe,[0 max(mean(spe,2))*3]);
% hold on
% scatter((c-500)/2,[1:1:size(spe,1)],5,'ow');hold off
colormap(hot)
axis([0 size(spe,2) 0.5 size(spe,1)+0.5])
set(gca,'xTick',[1:(size(spe,2)-1)/(size(fy,2)-1):size(spe,2)]);
set(gca,'xTicklabel',fy);
set(gca,'yTick',[]);

subplot('position',[0.46 0.07 0.5 0.2])
plot(mean(spe),'k')
colormap(hot)
axis([0 size(spe,2) 0 max(mean(spe))*2])
set(gca,'xTick',[1:(size(spe,2)-1)/(size(fy,2)-1):size(spe,2)]);
set(gca,'xTicklabel',fy);
set(gca,'yTick',[0:1000:round(max(mean(spe))*2/10)*10]);
ylabel('Intenisty')
xlabel('Wavelength (nm)')

%%

figure(2);
clf

subplot('position',[0.15 0.43 0.5 0.5]); hold on;
scatter(c,s,4,'bo','filled')
axis([600 800 0 1e4])
ylabel('Intenisty')
xlabel('Wavelength (nm)')
set(gca,'TickLength',[0.02,0.01]);
box on

subplot('position',[0.15 0.12 0.5 0.25])
[nb,xb]=hist(c,[600:1:800]); 
b=bar(xb,nb);
set(b,'facecolor',[1,1,1]*0.5);

set(gca,'TickLength',[0.02,0.01]);
xlabel('Wavelength (nm)')
ylabel('Count')
axis([600 800 0 1500])
std(c)

subplot('position',[0.7 0.43 0.25 0.5])
[nb,xb]=hist(s,[0:10:5e3]); 
nb=smooth(nb,3);
b=barh(xb,nb);
ylabel('Intensity ')
set(b,'facecolor',[1,1,1]*0.5);
set(gca,'TickLength',[0.02,0.01]);
xlabel('Count')
% set(gca,'YTickLabel','')
axis([0 2000 0 5000])

