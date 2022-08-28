%% read calibaration data
clear
clf
fnamein=sprintf('Tubulin+Alexa647_SR7nm 2018 April 26 18_34_47.spe');
readerobj = SpeReader(fnamein);
cali = squeeze(read(readerobj));
cali=mean(cali,3);
cali=mean(cali,1);
plot(cali)
%%
%Calibration data
xi=[487.7,546.5,611.6,707];
x0 = 26;
yi = [96,104,114,127];
res = (xi(end)-xi(1))/(yi(end)-yi(1));

%%
% pixel size: 160 nm
% Load spectrum file
clear
clc;
addpath('E:\Matlab_code\SPLM_code');
fname = 'Bao_AR1A1b_680+BRG1-647_cell7002';
fnamein=sprintf('%s.nd2',fname);
readerobj = SpeReader(fnamein);
A = read(readerobj);
A=double(squeeze(A));

% A=A-600;
%%
clear
clc;
addpath('E:\Matlab_code\SPLM_code');
fname = '6nm_Cell2002';
fnamein=sprintf('%s.nd2',fname);
A = SPLMload(fnamein,'nd2','double');
A = rot90(A,3);
%%
% Display

clf
for n=2:200
    imagesc(A(:,:,n),[300 2000]);
    title(n)
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    colormap(gray)
    set(gca,'color','red','Position',[0.05 0.1 0.9 0.8]);
    drawnow
%     pause(0.01)
end
%%
x1=11;
x2=84;
img=A(:,x1:x2,n);
imagesc(img(:,:,1));%axis equal;
%%
%Cropping image for ThunderStorm reconstruction
fnameout=sprintf('%s.dat',fname);
img=A(:,x1:x2,:);
fid1 = fopen(fnameout, 'w') ;
fwrite(fid1, img, 'uint16') ;
fclose(fid1);

imagesc(img(:,:,1));%axis equal;
%% Optional
%Calculate gobal background
% 
back=[];
for n=1:size(A,1)
    w=squeeze(A(n,220:320,:));
    w=mean(w,1);
    a=w<840;
    back(n,:)=mean(A(n,:,a),3);
end
imagesc(back)
%%
%Load reconstruction result from Thunderstorm
fnamein=sprintf('%s.csv',fname);
TR0=importdata(fnamein);
TR0=TR0.data;

TR=TR0(:,[3,4,2,6]);
PixelSize=160;
TR(:,1)=round(TR(:,1)/PixelSize+0.5);
TR(:,2)=round(TR(:,2)/PixelSize+0.5);

%Cropping for further process
w=3;
ax_min=TR(:,1)>w+3;
ax_max=TR(:,1)<max(TR(:,1))-w-1;
% ax_max=TR(:,1)<36;

ay_min=TR(:,2)>w+1;
ay_max=TR(:,2)<max(TR(:,2))-w-1;
% ay_max=TR(:,2)<36;
az_min=TR(:,3)>=w+1;
az_max=TR(:,3)<=max(TR(:,3))-w;
% az_max=TR(:,3)<=2000;
ai=TR(:,4)>=0;
a=ax_min.*ax_max.*ay_min.*ay_max.*az_min.*az_max.*ai;
a=logical(a);
TR=TR(a,:);
TR0=TR0(a,:);

out=hist3([TR0(:,3),TR0(:,4)],{0:10:max(TR0(:,3))*1.05,0:10:max(TR0(:,4)*1.05)});
imshow(out)

%%
clf;

for i=min(TR(:,3)):max(TR(:,3))
    temp=TR(TR(:,3)==i,:);
%     plot(temp(:,2),temp(:,1),'b+'); % x = temp(:,1) and y = temp(:,2)
    axis image;
    j=1;  
    hold on; 
    xy=temp;
    ind=[];
    while j<10
        hold on;
        peak = temp(find(xy(:,4)==max(xy(:,4))),1:2);
        xy(find(xy(:,4)==max(xy(:,4))),:) =00;
        plot(peak(2),peak(1),'r*'); % x = peak(1) and y = peak(2)
%         pause(1);
        ind(j,:)=peak(1);   
        j = j + 1;        
        
        if abs(xy(:,1) - peak(1)) < 2
            xy(xy(:,1)==peak(1),:) =0;  
        end
        if abs(ind(1:length(ind)-1) - peak(1))>1
            if abs(xy(:,1) - peak(1)) > 1
                    break;      
            end
        end        
    end
    
    if length(find(temp(:,1)==ind(end))) >1
        TRnew(i,:)= [0 0  0 0];
    elseif length(find(temp(:,1)==ind(end)))< 2
        TRnew(i,:)=temp(temp(:,1)==ind(end),:);   
    end
end


%%
TRnew(TRnew(:,1)==0,:)=[];
TR=TRnew;
%%
%optional; not recommend

%remove all overlapped events to collect reference spectra

FrameNumber=size(A,3);
Ta=[];
for n=1:FrameNumber
    w=TR(TR(:,3)==n,:);    
    a=ones(size(w,1),1);
    for m=1:size(w,1)
        tt=w(:,1); % y coordinates of the single molecules on same frame
        tt=abs (tt-tt(m));%
        tt(m)=150;% dimention of y
        at=tt>2; % difference between y coordinate
        a=a.*at;
    end
%     for mm=1:size(w,1)
%         tt=w(:,2);
%         tt=abs(tt-tt(m));
%         tt(mm)=50;
%         at2=tt>40;
%         a2=a.*at2;
%     end
    a=logical(a);
    Ta=[Ta;a];
end
Ta=logical(Ta);
TR=TR(Ta,:);
TR0=TR0(Ta,:);

%%

%Spectral calibration data (three wavelengths)
% slit = 43
xi=[487.7,546.5,611.6,707];
x0 = 47;
yi = [209,227,252,282];
fy = [600:25:800];%for xtick in plotting
fy2= [600:1:800];%for interpolation
coeffs=polyfit(xi, yi, 1);
fx = polyval(coeffs, fy);
fx2 = polyval(coeffs, fy2);
 
xc=x1-1;%cropping x-position
% xc=0;

y0=0;%Y shift

% input
s=TR(:,[2,1,3,4]);%this order relies on the output format of ThunderStorm
% s=TR(:,[2,1,3]);
% s=flipud(sortrows(s,4)); %sort by intensity
% bg=permute(bg,[2 1 3]);
s(:,1)=s(:,1)+xc;
s(:,2)=s(:,2);

aa=[];%Blinking image
bb=[];%Spectrum image
bb3=[];


as=[];%blinking profile
bs=[];%spectrum profile

bw=1;%blinking half width (actual width of the PSF=2*bw+1)

w=3;%half width of the calculated region
lbs=[-4:-3,3:4];%for local background subtraction

live=0;%live update: show=1; no show=0 

for n=1:size(s,1)


%     % %global background subtraction
    aa=A([-w:w]+s(n,2),[-w:w]+s(n,1),s(n,3))-back([-w:w]+s(n,2),[-w:w]+s(n,1));
    bb=A([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1),s(n,3))-back([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1));
%     aa2(:,:,n)=aa;  %
%     bb2(:,:,n)=bb;%
    as(n,:)=mean(aa(w+1-bw:w+1+bw,:));
    bs(n,:)=interp1([round(fx(1)):round(fx(end))],mean(bb(w+1-bw:w+1+bw,:)),fx2);
    bb3(:,:,n)=sum(bb(w+1-bw:w+1+bw,:));

%    % % local background subtraction
% sz=lbs+s(n,3);
% a=A([-w:w]+s(n,2),[-w:w]+s(n,1),s(n,3));
% aa=a-mean(A([-w:w]+s(n,2),[-w:w]+s(n,1),sz),3);
%   
% bb=b-mean(A([-w:w]+s(n,2)+y0,[round(fx(1)):round(fx(end))]-x0+s(n,1),sz),3);

%     as(n,:)=mean(aa(w+1-bw:w+1+bw,:));
%     bs(n,:)=interp1([round(fx(1)):round(fx(end))],mean(bb(w+1-bw:w+1+bw,:)),fx2);
%     bb3(:,:,n)=sum(bb(w+1-bw:w+1+bw,:));

    % show progress
    if mod(n*10000,size(s,1))<=50
        fprintf('Progress: %d \n',round(n/size(s,1)*100))
    end

    %Live display
    if live==1
        subplot('Position',[0.06 0.55 0.17 0.4])
            imagesc(aa,[0 500]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);

        subplot('Position',[0.3 0.55 0.65 0.4])
            imagesc(bb,[0 1000])
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
            axis([1 size(bs,2) -100 2000])
            set(gca,'XTick',[1:(size(bs,2)-1)/(size(fy,2)-1):size(bs,2)]);
            set(gca,'XTicklabel',fy);
        %     set(gca,'YTick',[]); 
        drawnow

    end
end
%%

%statistic analysis
% fy = [650:25:750];%for xtick in plotting
bs(isnan(bs))=0;
% bb3_1=squeeze(bb3)';
s=squeeze (sum(bb3,2));% after input from SPLM2, using bb3,2
% d=1e4;
% bs2=bs;
% bs2=bs(1+d:1e4+d,:);
% s=s(1+d:1e4+d)/100*3/0.925;
s=s/100*3/0.925;

x=500;
a1=s>x;
a2=s<=x+100000;
% a3=sum(bs(:,1:50),2)<=30000;
% a4=sum(bs(:,end-50:end),2)<=30000;
% a=logical(a1.*a2.*a3.*a4);
a=logical (a1.*a2);
s=s(a);
% spe=bs (a,:);
spe=bs(a,31:151);
% ns=na2(a,:);


%Calculate spectral centroid
wl=[630:1:750];
% wl=[500:2:800];
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
axis([0 125 0 max(mean(spe))*2])
set(gca,'xTick',[1:(size(spe,2)-1)/(size(fy,2)-1):size(spe,2)]);
set(gca,'xTicklabel',fy);
set(gca,'yTick',[0:1000:round(max(mean(spe))*2/10)*10]);
ylabel('Intenisty')
xlabel('Wavelength (nm)')

% a = c(c<=694.95+4&c>694.95-4);
% b = c(c<=705.16+4&c>705.16-4);
% d = c(c<=716.53+4&c>716.53-4);
% size(a),size(b),size(d)
std (c)
size (c)
a = c(c<=687+3&c>687-3);
b = c(c<=697+3&c>697-3);
% d = c(c<=710+3&c>710-3);
size(a)
size (b)
% mean (s),std(s)
%,size(b),size(d)


%% use gaussian fitting method instead of centroid method.
options = fitoptions('gauss1', 'Lower', [0 -inf 0 ],'Upper',[inf inf 30]);
for i=1:length(spe)
    warning off;
    f = fit(wl',spe(i,:)','gauss1',options);
    cx2(i,:) = f.b1;
end

% options = optimset('Display','off','TolFun',4e-16,'LargeScale','off');
% for i=1:length(spe)
% input = spe(i,:);
% data=input.*(input>0);
% cx1 = sum(wl.*(spe(i,:).^1))/sum(spe(i,:).^1);
% sx1 = sqrt(sum((abs((wl-cx1).^2).*data)/sum(data)));
% peak1 = max(max(data));
% 
% para = [cx1,sx1,peak1];
% fp1D = fminunc(@fitgaussian1D,para,options,data,x1);
% cx2(i) = fp1D(1);
% sx2 = fp1D(2);
% peak2 = fp1D(3);
% % wl2 = 0:0.01:wl(end);
% % output = abs(peak2)*(exp(-0.5*(wl2-cx2).^2./(sx2^2)));
% % plot(x2,output); xlim([500 800])
% end

c=cx2';
%%
figure(2);
clf

subplot('position',[0.15 0.43 0.5 0.5]); hold on;
scatter(c,s,4,'bo','filled')
axis([650 750 0 1e4])
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
axis([650 750 0 6000])
std(c)

% subplot('position',[0.7 0.43 0.25 0.5])
% [nb,xb]=hist(s,[0:10:2e3]); 
% nb=smooth(nb,3);
% b=barh(xb,nb);
% ylabel('Intensity ')
% set(b,'facecolor',[1,1,1]*0.5);
% set(gca,'TickLength',[0.02,0.01]);
% xlabel('Count')
% set(gca,'YTickLabel','')
% axis([0 1000 0 2000])
%% calcaulte FWHM
cs = [round(c)' round(s)];
cs = sortrows(cs,1);
cs2 = cs;

for i = 1:length(cs2)
   cs2(cs2(:,1)==cs2(i,1),2)=max(cs2(cs2(:,1)==cs2(i,1),2));
end
% cs2(cs2(:,1)>(f.b1+3*f.c1),:)=[]; % remove out of 3 sigma value; this is only for calculating fwhm and better visualization.
base = min(cs2(:,2));

subplot('position',[0.15 0.43 0.5 0.5]); hold on;
scatter(cs2(:,1),cs2(:,2),4,'ko','filled')
axis([500 800 0 7e3])
ylabel('Intenisty')
xlabel('Wavelength (nm)')
set(gca,'XTickLabel','')
set(gca,'TickLength',[0.02,0.01]);
box on

f = fit(cs2(:,1),cs2(:,2)-base,'gauss1');
yfit =  f.a1.*exp(-((cs2(:,1)-f.b1)./f.c1).^2) + base;
fwhm = f.c1*2.35

% subplot('position',[0.15 0.43 0.5 0.5]); hold on;
plot(cs2(:,1),yfit,'g-','Linewidth',1);
%%
% %output

bb3=squeeze (bb3);
save1 = sprintf('%s_TR',fname);
save2 = sprintf('%s_TR0',fname);
save3 = sprintf('%s_spe',fname);
save4 = sprintf('%s_bb3',fname);

save(save1,'TR','-ascii') %position in formation
save(save2,'TR0','-ascii') %position information
save(save3,'bs','-ascii') %spectral information
save(save4,'bb3','-ascii') %spectral information,raw


%%
Data_new=A(1:end,200:230,:);
mu = mean(Data_new,3);
for i=1:500
    Data_new(:,:,i)=Data_new(:,:,i)-mu; 
end

Count_Readout1 = std(Data_new(:),1)  % Offset 패턴 노이즈가 제거 + 에버리지: readout noise in counts 

% method 2
data1 = Data_new(:,:,1);
data2 = Data_new(:,:,2);
data = data1 - data2;
Count_Readout2 = std(data(:))/sqrt(2)

 Nr = Count_Readout1*3/100
 
 %%
 tic
 options = optimset('Display','off','TolFun',4e-16,'LargeScale','off');
for n = 1:size(s,1)
    m = aa2(:,:,n);
    [sizey sizex] = size(m);
    [cx,cy,sx,sy] = fit_Center(m);
    amp = max(m(:));
    mx2 = sum(m);
    x1D = 1:sizex;
    ip1D = [cx,sx,amp];
    try
        fp1D = fminunc(@fitgaussian1D,ip1D,options,mx2,x1D);
    end
        cx(:,n) = fp1D(1);
        sx(:,n) = fp1D(2);
        amp2 = fp1D(3);
        x = 1:sizex;
        mfit = amp2*(exp(-0.5*(x-cx(:,n)).^2./(sx(:,n)^2)));
        aa3(n,:) = mfit;
end
toc;
aa3(sx>10,:)=[];
aa3(cx<0,:)=[];
aa4=mean(aa3);
aa4=aa4/max(aa4);

%%
subplot 211
n=10;
 plot(bs(n,:))
subplot 212
bs2=fdeconv(bs(n,:),aa4)
plot(bs2);
%  aa5=conv(bs2,aa4);
%  plot(aa5)
 

%%
 tic;
for j=1:1000
    try        
        [cx(j), cy(j), sx(j), sy(j), Ip(j)] = fit_2DGaussian(aa2(:,:,j));    
        [sizey sizex] = size(aa2(:,:,j));
        [x,y] = meshgrid(1:sizex,1:sizey);
        m_fit(:,:,j) = abs(Ip(j))*(exp(-0.5*(x-cx(j)).^2./(sx(j)^2)-0.5*(y-cy(j)).^2./(sy(j)^2)));
    end
end
m_fit=mean(m_fit,3);
m_fit=m_fit/max(m_fit(:));
subplot(1,2,1); imagesc(aa2(:,:,j));axis image; set(gca,'YTick',[]);  set(gca,'XTick',[]);
subplot(1,2,2); imagesc(m_fit);  colormap(gray); axis image; set(gca,'YTick',[]);  set(gca,'XTick',[]);
toc;

%%
subplot 121
imagesc(bb2(:,:,1)); plot(sum(bb2(:,:,1)))
% subplot 122
hold on;
bb3=deconvblind(bb2,m_fit);
imagesc(bb3(:,:,1)); plot(sum(bb3(:,:,1)))

for n=1:size(s,1)
bs(n,:)=interp1([round(fx(1)):round(fx(end))],mean(bb3(w+1-bw:w+1+bw,:,n)),fx2);
end