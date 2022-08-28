%%
% Load spectrum file
clear
clc;
fnamein=sprintf('Alexa546_Tubulin_wide_2018 March 22 11_29_13.spe');

readerobj = SpeReader(fnamein);

A = read(readerobj);

A=double(squeeze(A));
%%
% Display

clf
for n=1:500
    imagesc(A(:,:,n),[500 3000]);
    title(n)
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    colormap(gray)
    set(gca,'color','red','Position',[0.05 0.1 0.9 0.8]);
    drawnow
%     pause(0.01)
end

%%
%Cropping image for ThunderStorm reconstruction

fnameout=sprintf('Alexa546_Tubulin_wide_2018 March 22 11_29_13.dat');
img=A(:,37:113,:);
fid1 = fopen(fnameout, 'w') ;
fwrite(fid1, img, 'uint16') ;
fclose(fid1);

%%
%Optional
%Calculate gobal background

back=[];
for n=1:size(A,1)
    w=squeeze(A(n,20:78,:));
    w=mean(w,1);
    a=w<1185;
    back(n,:)=mean(A(n,:,a),3);
end

imagesc(back)
%%

%Load reconstruction result from Thunderstorm
TR0=importdata('Alexa546_Tubulin_wide_2018 March 22 11_29_13.csv');
TR0=TR0.data;

TR=TR0(:,[2,3,1,5]); % x,y,frame number and intensity
% TR=flipud(round(sortrows(TR,4)));

PixelSize=160;
TR(:,1)=round(TR(:,1)/PixelSize+0.5);
TR(:,2)=round(TR(:,2)/PixelSize+0.5);

%Cropping for further process
w=3;
ax_min=TR(:,1)>1+w;
ax_max=TR(:,1)<120-w-1;
ay_min=TR(:,2)>w+1;
ay_max=TR(:,2)<77-w-1;
az_min=TR(:,3)>=10;
az_max=TR(:,3)<=5000-10;
ai=TR(:,4)>=0;
a=ax_min.*ax_max.*ay_min.*ay_max.*az_min.*az_max.*ai;
a=logical(a);
TR=TR(a,:);
TR0=TR0(a,:);

out=hist3([TR0(:,2),TR0(:,3)],{0:10:12000,0:10:20000});
imshow(out)

%%
%optional

%remove all overlapped events to collect reference spectra

FrameNumber=size(A,3);
Ta=[];
for n=1:FrameNumber
    w=TR(TR(:,3)==n,:);    
    a=ones(size(w,1),1);
    for m=1:size(w,1)
        tt=w(:,1);
        tt=tt-tt(m);
        tt(m)=100;
        at=tt>1;
        a=a.*at;
    end
    a=logical(a);
    Ta=[Ta;a];
end
Ta=logical(Ta);
TR=TR(Ta,:);
TR0=TR0(Ta,:);

%%

%Spectral calibration data (three wavelengths)
yi=[227,245,266,297];
xi=[487.7,546.5,611.6,707];


fy = [525:25:675];
fy2= [500:1:800];%for interpolation
coeffs=polyfit(xi, yi, 1);
fx = polyval(coeffs, fy);
fx2 = polyval(coeffs, fy2);
%%

x0=73;%slit position
xc=36;%cropping x-position

y0=2;%Y shift

% input
s=TR(:,[2,1,3,4]);%this order relies on the output format of ThunderStorm
% s=flipud(sortrows(s,4)); %sort by intensity

s(:,1)=s(:,1)+xc;
s(:,2)=s(:,2);

aa=[];%Blinking image
bb=[];%Spectrum image

as=[];%blinking profile
bs=[];%spectrum profile

bw=1; %blinking half width (actual width of the PSF=2*bw+1)

w=3;%half width of the calculated region
lbs=[-6:-3,3:6];%for local background subtraction

live=0;%live update: show=1; no show=0 

for n=1:size(s,1)
    
% % % local background subtraction
% sz=lbs+s(n,3);
% a=A([-w:w]+s(n,2),[-w:w]+s(n,1),s(n,3));
% aa=a-mean(A([-w:w]+s(n,2),[-w:w]+s(n,1),sz),3);
% b=A([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1),s(n,3));
% bb=b-mean(A([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1),sz),3);

%global background subtraction
aa=A([-w:w]+s(n,2),[-w:w]+s(n,1),s(n,3))...
    -back([-w:w]+s(n,2),[-w:w]+s(n,1));
bb=A([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1),s(n,3))...
    -back([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1));
% 
as(n,:)=mean(aa(w+1-bw:w+1+bw,:));
bs(n,:)=interp1([round(fx(1)):round(fx(end))],mean(bb(w+1-bw:w+1+bw,:)),fx2);

% show progress
if mod(n*10000,size(s,1))<=50
    fprintf('Progress: %d \n',round(n/size(s,1)*100))
end

%Live display
    if live==1
        subplot('Position',[0.06 0.55 0.17 0.4])
            imagesc(aa,[0 4000]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);

        subplot('Position',[0.3 0.55 0.65 0.4])
            imagesc(bb,[0 3000])
            set(gca,'XTick',[1:(size(bb,2)-1)/(size(fy,2)-1):size(bb,2)]);
            set(gca,'XTicklabel',[540:20:700]);
            set(gca,'YTick',[]);
            colormap(gray)

        subplot('Position',[0.06 0.1 0.17 0.35])
            plot(as(n,:),'k');
            % hold on
            % plot(smooth(bs(n,:),20),'r');hold off
            axis([1 size(as,2) -50 4000])
            set(gca,'XTick',[1:(size(bs,2)-1)/(size(fy,2)-1):size(bs,2)]);
            set(gca,'XTicklabel',fy);
        %     set(gca,'YTick',[]); 

        subplot('Position',[0.3 0.1 0.65 0.35])
            plot(bs(n,:),'k');
            % hold on
            % plot(smooth(bs(n,:),20),'r');hold off
            axis([1 301 -100 4000])
%             set(gca,'XTick',[1:(size(bs,2)-1)/(size(fy,2)-1):size(bs,2)])
            set(gca,'XTick',[1:50:301]);
            set(gca,'XTicklabel',[500:50:800]);
        %     set(gca,'YTick',[]); 
        drawnow
    end

end
%%
% %output
% 
save('ATPB+Tub#5#secondcell 21_11_44_TR','TR0','-ascii') %position information
save('ATPB+Tub#5#secondcell 21_11_44_spe','bs','-ascii') %spectral information
% 
%%

%load

%Calibration data
yi=[208,225,248];
xi=[487.7,546.5,611.6];
fy = [550:25:750];%for xtick in plot
fy2= [500:2:750];%for interpolation
coeffs=polyfit(xi, yi, 1);
fx = polyval(coeffs, fy);  
fx2 = polyval(coeffs, fy2);
%%

bs=load('ATPB+Tub#5#secondcell 21_11_44_spe');
na2=load('ATPB+Tub#5#secondcell 21_11_44_TR');

%%

%statistic analysis

bs(isnan(bs))=0;
s=sum(bs,2);
s=s/100*3/0.95;

%Flitering based on spectral intensity
a1=s>=2000;
a2=s<=1e5;
% a3=sum(bs(:,1:80),2)<=1e3;
% a4=sum(bs(:,end-30:end),2)<=2e5;
% a=logical(a1.*a2.*a3.*a4);
a=logical (a1.*a2);
s=s(a);
spe=bs(a,:);
% ns=na2(a,:);


%dispersion correction

%replace corr with dispersion calibration data
% corr=meshgrid(linspace(0.7,0.55,451),linspace(1,1,size(spe,1)));
% spe=spe./corr;
%%
clf
wavel=500:700;
cfit=zeros (1,size(spe,1));
r2=zeros (1,size(spe,1));
for n = 1:size (spe,1)
    [ft,gn] = fit(wavel',spe(n,1:201)','gauss1');
%     plot (ft);
% %     axis ([500 700 -100 1500]);
%     title (n);
%     hold on
%     plot (wavel,spe(n,1:201));
%     pause (0.1)
%     axis ([500 700 -100 1500]);
%     hold off
    cfit (n)= ft.b1;
    r2 (n)=gn.rsquare;
end
    cfit1=cfit(r2>=0.5);
    sfit=s(r2>=0.5);
%%
clf

 
%Calculate spectral centroid
wl=[500:1:800];
c=[];
for n=1:size(spe,1)
    c(n)=sum(wl.*(spe(n,:).^1))/sum(spe(n,:).^1);
end 

subplot('position',[0.12 0.35 0.28 0.6])

barh(flipud(sum(spe,2)),'EdgeColor','None')
axis([0 1e5 0.5 size(spe,1)+0.5])
ylabel('Blinking')
xlabel('Intenisty')

subplot('position',[0.46 0.35 0.5 0.6])
imagesc(spe (:,25:175),[500 2500]);
% hold on
% scatter((c-500)/2,[1:1:size(spe,1)],5,'ow');hold off
colormap(hot)
axis([0 150 0.5 size(spe,1)+0.5])
set(gca,'xTick',[0:25:150]);
set(gca,'xTicklabel',[525:25:675]);
set(gca,'yTick',[]);

subplot('position',[0.46 0.07 0.5 0.2])
plot(mean(spe),'k')
colormap(hot)
axis([0 size(spe,2) 0 2000])

set(gca,'xTick',[1:50:301]);
set(gca,'xTicklabel',[500:50:800]);
set(gca,'yTick',[0:400:4000]);
ylabel('Intenisty')
xlabel('Wavelength (nm)')

%%
clf

subplot('position',[0.15 0.43 0.5 0.5])

% scatplot(c,s,'voronoi',5,100,1,1,2);
% colormap(jet)
scatter(c,s,4,'ko','filled')
axis([550 650 800 5e3])
ylabel('Intenisty')
xlabel('Wavelength (nm)')
set(gca,'XTickLabel','')
set(gca,'TickLength',[0.02,0.01]);
box on


subplot('position',[0.15 0.12 0.5 0.25])
[nb,xb]=hist(c,[500:5:700]); 
b=bar(xb,nb);
set(b,'facecolor',[1,1,1]*0.5);

set(gca,'TickLength',[0.02,0.01]);
xlabel('Wavelength (nm)')
ylabel('Count')
axis([550 650 0 1000])
% std(c)

subplot('position',[0.7 0.43 0.25 0.5])
[nb,xb]=hist(s,[0:1e2:2e3]); 
nb=smooth(nb,3);
b=barh(xb,nb);
ylabel('Intensity ')
set(b,'facecolor',[1,1,1]*0.5);
set(gca,'TickLength',[0.02,0.01]);
xlabel('Count')
set(gca,'YTickLabel','')
axis([0 5000 0 2e3])

% result = [c',s];
% xlswrite ('RB 2017 August 31 16_19_00.xls',result,1,'A1');