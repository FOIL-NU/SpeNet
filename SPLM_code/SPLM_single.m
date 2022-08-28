% Load spectrum file
clear
clc;
fnamein=sprintf('RhB_PS_0.1nM__time_50ms 2017 December 09 20_00_09.spe');

readerobj = SpeReader(fnamein);

A = read(readerobj);

A=double(squeeze(A));


%%
% draw an image
clf
for n=1:500
    imagesc(A(:,:,n),[400 2000]);
    title(n)
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    colormap(gray)
    set(gca,'color','none','Position',[0.05 0.1 0.9 0.8]);
    drawnow
    pause(0.01)
end

%% 

back=[];
for n=1:size(A,1)
   w=squeeze(A(n,26:65,:));
   w=mean(w,1);
   a=w<710; 
   back(n,:)=mean(A(n,:,a),3);
end

imagesc(back)
%%
back=repmat(back,[1,1,5000]);
A=A-back;
A=A(6:120,:,:);
%%
%Spectral calibration data (three wavelengths)
yi=[211,229,252,286]+3;
xi=[487.7,546.5,611.6,707];

fy = [525:25:675];%for xtick in plotting
fy2= [525:1:675];%for interpolation
coeffs=polyfit(xi, yi, 1);
fx = polyval(coeffs, fy);
fx2 = polyval(coeffs, fy2);

x0=47;% pixel number determining the localization postion . = slit position
xc=25;%cropping x-position

y0=-1;%Y shift

% input
s=[45,48]; %s=[31,17];%position
sz=[120:200]';
s=[repmat(s,size(sz)),sz];

aa=[];%Blinking image
bb=[];%Spectrum image

as=[];%blinking profile
bs=[];%spectrum profile

bw=1;%blinking half width (actual width of the PSF=2*bw+1)

w=3;%half width of the calculated region
lbs=[-6:-3,3:6];%for local background subtraction

live=1;%live update: show=1; no show=0 

for n=1:size(s,1)
    
%local background subtraction
% sz=lbs+s(n,3);
% a=A([-w:w]+s(n,1),[-w:w]+s(n,2),s(n,3));
% aa=a-mean(A([-w:w]+s(n,1),[-w:w]+s(n,2),sz),3);
% b=A([-w:w]+s(n,1)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,2),s(n,3));
% bb=b-mean(A([-w:w]+s(n,1)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,2),sz),3);

% %global background subtraction
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

%Flitering based on spectral intensity
a1=s>=60000;
a2=s<=1e7;
% a3=sum(bs(:,1:90),2)<=5e3;
% a4=sum(bs(:,end-25:end),2)<=5e3;
a=logical(a1.*a2);
s=s(a);
spe=bs(a,:);
% ns=na2(a,:);


%dispersion correction

%replace corr with dispersion calibration data
% corr=meshgrid(linspace(0.7,0.55,451),linspace(1,1,size(spe,1)));
% spe=spe./corr;


%Calculate spectral centroid
wl=[525:1:675];
c=[];
for n=1:size(spe,1)
    c(n)=sum(wl.*(spe(n,:).^1))/sum(spe(n,:).^1);
end 

subplot('position',[0.12 0.35 0.28 0.6])

barh(flipud(sum(spe,2)),'EdgeColor','None'  )
axis([1e4 1e5 0.5 size(spe,1)+0.5])
ylabel('Blinking')
xlabel('Intenisty')

subplot('position',[0.46 0.35 0.5 0.6])
imagesc(spe, [100 1600]);
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
axis([0 size(spe,2) 0 800])

set(gca,'xTick',[1:(size(spe,2)-1)/(size(fy,2)-1):size(spe,2)]);
set(gca,'xTicklabel',fy);
set(gca,'yTick',[0:400:4000]);
ylabel('Intenisty')
xlabel('Wavelength (nm)')

%%
%%output to excel  
result = [spe];
xlswrite ('PMMA_4100_5000_24_76.xls',result,1,'A1');

%%
clf

subplot('position',[0.15 0.43 0.5 0.5])

% scatplot(c,s,'voronoi',5,100,1,1,2);
% colormap(jet)
scatter(c,s,4,'ko','filled')
axis([560 600 2e4 7e4])
ylabel('Intenisty')
% xlabel('Wavelength (nm)')
set(gca,'XTickLabel','')
set(gca,'TickLength',[0.02,0.01]);


box on


subplot('position',[0.15 0.12 0.5 0.25])
[nb,xb]=hist(c,[400:1:800]); 
b=bar(xb,nb);
set(b,'facecolor',[1,1,1]*0.5);

set(gca,'TickLength',[0.02,0.01]);
xlabel('Wavelength (nm)')
ylabel('Count')
axis([560 600 0 20])
sd=sprintf('S.D.=%.2f nm',std(c));
text(580,16,sd)

subplot('position',[0.7 0.43 0.25 0.5])
[nb,xb]=hist(s,[0:0.2e4:6e5]); 
nb=smooth(nb,3);
b=barh(xb,nb);
set(b,'facecolor',[1,1,1]*0.5);
set(gca,'TickLength',[0.02,0.01]);
xlabel('Count')
set(gca,'YTickLabel','')
axis([0 20 2e4 7e4])

