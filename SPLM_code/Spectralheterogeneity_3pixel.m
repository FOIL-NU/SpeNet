%%
addpath ('E:\Matlab_code\SPLM_code'); 
% Load spectrum file
clear
clc;
fnamein=sprintf('RhB_0.1nM_Glass 2017 December 08 16_56_30.spe');
readerobj = SpeReader(fnamein);

A = read(readerobj);

A=double(squeeze(A));
%%
A=A (:,:,2000:end);
%%
clf
for n=1:500
    imagesc(A(:,:,n),[600 1500]);
    title(n)
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    colormap(gray)
    set(gca,'color','none','Position',[0.05 0.1 0.9 0.8]);
    drawnow
    pause(0.01)
end

%%
%Calculate gobal background

back=[];  
for n=1:size(A,1)
    w=squeeze(A(n,202:270,1:end)); %why 242 to 270?
    w=mean(w,1);
    a=w<610;
    back(n,:)=mean(A(n,:,a),3);
end
clf
imagesc(back)
%%
% back=mean(A(:,:,100:150),3);
% imagesc(back)
% backall = load ('back.mat','-mat');
% back = backall.back;
back=repmat(back,[1,1,size(A,3)]);
A=A-back;
A=A(10:120,:,:);
%%
yi=[211,229,252,286];
xi=[487.7,546.5,611.6,707];
slit = 46;
fy = [540:25:690];
fy2= [500:1:800];
coeffs=polyfit(xi, yi, 1);
fx = polyval(coeffs, fy);
fx2 = polyval(coeffs, fy2);

% slit=45;

bb=[];


% s=[38,103,2];
as=[];
bs=[];
rs=[];
%%
back1 = mean(back,3);
x_loc =[];
for n=1:120
    [x,xlocation] = max (back(n,1:100));
    x_loc = [x_loc,xlocation];
end
x_loc= x_loc (x_loc>slit-10&x_loc<slit+10);
slit_shift = max (x_loc)-min(x_loc);   

%%
as = [];
bs = [];
s = [];
bsi = [];
fx1 = round (fx(1:end));
for n=1:size (A,3);
    frame = A (:,:,n);
    as = frame (:,slit-20:slit+20);
    s=sum(as,2);
    t=1000;
    zeroders (:,n) = s';
    [~,locs_Rwave] = findpeaks (s,'MinPeakHeight',t,'MinPeakDistance',5);
    if size (locs_Rwave,2)>=1
        for m=1:size(locs_Rwave,1)
            [~,I] = max (frame(locs_Rwave(m),1:100));
            d=I-slit;
            single_spetrum_image = frame (locs_Rwave(m)-1:locs_Rwave(m)+1,fx1(1)+d:fx1(end)+d);
            bs = [bs;sum(single_spetrum_image,1)];
            bsi = [bsi,sum(sum(single_spetrum_image))];
        end
    end
end
%%
bs3=[];
for n=1:size(bs,1)
    bs3(n,:)=interp1([round(fx(1)):round(fx(end))],bs(n,:),fx2);
end

bs3(isnan(bs3))=0;
%%
% bs3(isnan(bs3))=0;
s=bsi;
s = s/100*3/0.95*1.5;
%Flitering based on spectral intensity
a1=s>1000;
a2=s<=1e6;
% a3=sum(bs3(:,1:20),2)<=3e3;
% a4=sum(bs3(:,end-10:end),2)<=2e3;
% a=logical(a1.*a2.*a3.*a4);
a=logical(a1.*a2);
s=s(a);
spe=bs3(a,:);
c = [];
for n = 1:size (spe,1)
    wavelength = 500:1:800;
    c1 = sum (wavelength.*spe (n,:))/sum (spe(n,:));
    c = [c,c1];
end
   
%%
clf
% s=s/300*4.2;

% a1=s>=1e1;
% a2=s<=1e5;
% % a3=sum(bs(:,1:50),2)<=5e2;
% % a4=sum(bs(:,end-80:end),2)<=10e2;
% a=logical(a1.*a2);
% spe2=spe(a,:);
% c2=c(a);

subplot('position',[0.1 0.15 0.25 0.81])
imagesc(spe,[0 2000])
% axis([500 700 0 200])
ylabel('Photons')
xlabel('Wavelength (nm)')
% set(gca,'XTickLabel','')
set(gca,'TickLength',[0.02,0.01]);
set(gca,'XTick',[600:100:800]);
set(gca,'XTicklabel',[600:100:800]);
box on

subplot('position',[0.45 0.51 0.3 0.45])
scatter(c,s,3,'ko','filled'); hold on
% plot([555,625],[1,1]*0.25e4,'r');
% plot([555,625],[1,1]*0.3e4,'r');
% plot([555,555],[0.25,0.3]*1e4,'r');
% plot([605,625],[0.25,0.3]*1e4,'r');
hold off
axis([600 800 0 1e4]);
ylabel('Photons');
xlabel('Wavelength (nm)');
% set(gca,'XTickLabel','')
set(gca,'TickLength',[0.02,0.01]);
set(gca,'xTick',[600:25:800])
box on

colormap(hot)

subplot('position',[0.45 0.15 0.3 0.25])
[nb,xb]=hist(c,[600:3:800]); 
b=bar(xb,nb/sum(nb)*100);
set(b,'facecolor',[1,1,1]*0.5);

set(gca,'TickLength',[0.02,0.01]);
xlabel('Wavelength (nm)')
ylabel('Probability %')
axis([600 800 0 40])
set(gca,'xTick',[600:25:800])

mean(s)
subplot('position',[0.8 0.51 0.15 0.45])
[nb,xb]=hist(s,[0:100:10000]); 
% nb=smooth(nb,3);
b=barh(xb,nb/sum(nb)*100);
set(b,'facecolor',[1,1,1]*0.5);
set(gca,'TickLength',[0.02,0.01]);
xlabel('Probability %')
set(gca,'YTickLabel','')
axis([0 20 0 5000])

% %%output to excel  
% result = [c',s];
% xlswrite ('PMMA 2017 September 28 16_29_41.xls',result,1,'A1');