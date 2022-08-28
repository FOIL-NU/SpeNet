%%
addpath ('E:\Matlab_code\SPLM_code'); 
% Load spectrum file
clear
clc;
fnamein=sprintf('RhB_only 2017 December 10 15_04_15.spe');
readerobj = SpeReader(fnamein);

A = read(readerobj);

A=double(squeeze(A));
%%
A=A(:,:,1:2000);
%%
Ave1= mean (A ,3);
imagesc (Ave1);
slit_position = [];
for n = 1: size (Ave1,1)
    [m,I] = max (Ave1(n,1:100));
    slit_position = [slit_position, I];
end
slit_position= slit_position';
%%
clf
for n=2000:2220
    imagesc(A(:,:,n),[400 2000]);
    title(n)
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    colormap(gray)
    set(gca,'color','none','Position',[0.05 0.1 0.9 0.8]);
    drawnow
%     pause(0.01)
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
A=A(50:90,:,:);
%%
yi=[210,229,252,286]+3;
xi=[487.7,546.5,611.6,707];
slit=49
fy = [525:25:675];
fy2= [500:1:800];
coeffs=polyfit(xi, yi, 1);
fx = polyval(coeffs, fy);
fx2 = polyval(coeffs, fy2);

x0=46;

bb=[];


% s=[38,103,2];
as=[];
bs=[];
rs=[];

%%
for n=1:size(A,3)

b=A(:,round(fx(1)):round(fx(end)),n);
s=sum(b,2);
t=4000;
% a=s>=t;
% 
as(:,n)=s';

warning('off')
[~,locs_Rwave] = findpeaks(s,'MinPeakHeight',t,'MinPeakDistance',5);
if size(locs_Rwave,2)>=1
    for m=1:size(locs_Rwave,1)
        
        bs=[bs;sum(b(locs_Rwave(m)-1:locs_Rwave(m)+1,:),1)];
        rs=[rs;[n,locs_Rwave(m)]];
    end
   
end

%     subplot 211
% 
%     plot(s);hold on
%     plot(locs_Rwave,s(locs_Rwave));hold off
% 
%     %     imagesc(b,[0 1000])
%     %     set(gca,'XTick',[]);
%     % %     set(gca,'XTicklabel',fy);
%     %     set(gca,'YTick',[]);   
%     % %     set(gca,'color','none','Position',[0.3 0.55 0.65 0.4]);
%     % colormap(gray)
%     % 
%     subplot 212
% 
%         imagesc(bs,[0 1000])
%         set(gca,'XTick',[]);
%     %     set(gca,'XTicklabel',fy);
%         set(gca,'YTi  ck',[]);   
%     %     set(gca,'color','none','Position',[0.3 0.55 0.65 0.4]);
%     colormap(gray)

end

clf
imagesc(as,[000 10000]);hold on
plot(rs(:,1),rs(:,2),'wo');hold off
l=[1];
d=200;
rs2=rs;
bs2=bs;  
% for m=1:size(l,2)
%     a1=rs2(:,1)>(l(m)-d);
%     a2=rs2(:,1)<(l(m)+d);
%     a=logical(a1.*a2);
%     rs2(a,:)=[];
%     bs2(a,:)=[];
% end

%%
clf

bs2(isnan(bs2))=0;

bs3=[];
for n=1:size(bs2,1)
    bs3(n,:)=interp1([round(fx(1)):round(fx(end))],bs2(n,:),fx2);
end

bs3(isnan(bs3))=0;

% save('QR_spe','bs2','-ascii')

sig=sum(bs3,2);


a1=sig>1e0;
a2=sig<1e10;
a=a1.*a2;

bs3=bs3(logical(a),:);
bs3=bs3(randperm(size(bs3,1)),:);
% subplot 212

imagesc(bs3,[00 6000])
colormap(hot)
axis([0 200 0.5 size(bs3,1)+0.5  ])

ylabel('Blinking')
xlabel('Wavelength (nm)')
%%
save('NR_0.6nM 2018 February 04 18_27_02','bs','-ascii')

%%
bs1=load('NR_0.6nM 2018 February 04 18_27_02');
bs2=load('NR_0.6nM 2018 February 04 18_28_58');
% bs3=load('ATTO655_0.6nM_Glass 2018 February 05 14_16_56');
% bs1=bs3.bs3;
%%
bs=cat (1,bs1,bs2);

%%
% bs3(isnan(bs3))=0;
s=sum(bs,2);
s = s/100*3/0.95;
a1=s>=200;
a2=s<=1e5;
% a3=sum(bs3(:,1:20),2)<=3e3;
% a4=sum(bs3(:,end-10:end),2)<=2e3;
% a=logical(a1.*a2.*a3.*a4);
a=logical(a1.*a2);
s=s(a);
spe=bs3(a,:);


wl=[500:1:800];
c=[];
for n=1:size(spe,1)
    c(n)=sum(wl.*spe(n,:))/sum(spe(n,:));
end
% w=s./max(spe')';
% 
% a1=w>=0;
% a2=w<=31;
% a=logical(a1.*a2);
% s=s(a);
% c=c(a);
% spe=spe(a,:);

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
imagesc(spe,[200 1500])
% axis([500 700 0 200])
ylabel('Photons')
xlabel('Wavelength (nm)')
% set(gca,'XTickLabel','')
set(gca,'TickLength',[0.02,0.01]);
set(gca,'XTick',[500:100:700]);
set(gca,'XTicklabel',[500:100:700]);
box on

subplot('position',[0.45 0.51 0.3 0.45])
scatter(c,s,3,'ko','filled'); hold on
% plot([555,625],[1,1]*0.25e4,'r');
% plot([555,625],[1,1]*0.3e4,'r');
% plot([555,555],[0.25,0.3]*1e4,'r');
% plot([605,625],[0.25,0.3]*1e4,'r');
hold off
axis([550 650 0 5e3])
ylabel('Photons')
xlabel('Wavelength (nm)')
% set(gca,'XTickLabel','')
set(gca,'TickLength',[0.02,0.01]);
set(gca,'xTick',[500:25:700])
box on

colormap(hot)

subplot('position',[0.45 0.15 0.3 0.25])
[nb,xb]=hist(c,[500:5:800]); 
b=bar(xb,nb/sum(nb)*100);
set(b,'facecolor',[1,1,1]*0.5);

set(gca,'TickLength',[0.02,0.01]);
xlabel('Wavelength (nm)')
ylabel('Probability %')
axis([550 650 0 40])
set(gca,'xTick',[525:25:675])

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


%%
clf
[nb,xb]=hist(c,[500:2.5:800]); 
% x=x';y=y';aaa = mean (spe)';
% plot (x,y)
b = bar (xb, nb/sum(nb)*100,1);
set(b,'facecolor',[0.7,0.7,0.7]);
axis ([525 675 0 40]);
txt1 =[ 's = ',num2str(round (std(c),1)),' nm'];
txt2 = [fnamein];
text (640, 35, txt1,'FontSize',15);
text (540, 35, txt2,'FontSize',15);
set (gca,'xTick', ([525:50:675]));
set (gca,'yTick', ([0:10:40]));
set (gca, 'FontSize',15);
box on
% xlabel('Wavelength (nm)');
% ylabel('Probability %');
%%
ft= fit (xb',(nb/sum(nb)*100)','gauss1')'
hold on
plot (ft);
%%
clf
subplot('position',[0.2 0.46 0.7 0.5])
n=10;
sp=spe(n,:);
[fitobject,gof,output] = fit([1:301].',sp.','gauss1');
sp=sp'+output.residuals;

plot(sp)

for n=1:100
c2(n)=randn*0.5;
sp2(:,n)=circshift(sp,round(c2(n)))+(sp+5000).*randn(301,1)/60;
end
imagesc(sp2',[0 6500]);hold on
scatter(flipud(c2)+124,[1:1:200],2,'ko');hold off

axis([75 175 0 100])
ylabel('Photons')
xlabel('Wavelength (nm)')
% set(gca,'XTickLabel','')
set(gca,'TickLength',[0.02,0.01]);
set(gca,'XTick',[75:25:175]);
set(gca,'XTicklabel',[525:25:625]);
% box on

subplot('position',[0.2 0.1 0.7 0.2])
[nb,xb]=hist(c2+574,[525:0.5:625]);

std(c2)
% nb=smooth(nb,3);
b=bar(xb,nb/sum(nb)*100);
set(b,'facecolor',[1,1,1]*0.5);
set(gca,'TickLength',[0.02,0.01]);
ylabel('Probability %')
% set(gca,'YTickLabel','')
set(gca,'XTick',[550:5:600]);
axis([565 580 0 40])

%% find emission max

cc = [];
for n=1:size(spe,1)
%     max ()
y = spe (n,:);
xIndex = find(y == max(y), 1, 'first');
maxXValue = wl(xIndex);
%     plot (wl,spe(n,:));
    width = 20;
   wavel= maxXValue-width:maxXValue+width;
      ft =fit (wavel',(spe(n,wavel(1:end)-500)'),'gauss1');
%            hold on
%       plot (ft);
      cc = [cc,ft.b1];
%     drawnow
%     pause (0.1);
%     hold off
end
%%
%scatter (cc1,s,3,'filled','k');
a1 = cc>500;
a2 = cc<800;
a=logical (a1.*a2);
cc1=cc(a);
s1=s(a);
[nb,xb] = hist (cc1,[500:4.5:800]);
bar (xb,nb/sum(nb)*100);
axis ([525 675 0 30]);
ft = fit (xb',(nb/sum(nb)*100)','gauss1');
hold on
plot (ft);

