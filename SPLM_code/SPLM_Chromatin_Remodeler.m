%%
clear
clc;
addpath('E:\Matlab_code\SPLM_code');
fname = 'ell1_008';
fnamein=sprintf('%s.nd2',fname);
A = SPLMload(fnamein,'nd2','double');
A = rot90(A,3);

%%
% Display

clf
for n=300:500
    imagesc(A(:,:,n),[400 2500]);
    title(n)
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    colormap(gray)
    set(gca,'color','red','Position',[0.05 0.1 0.9 0.8]);
    drawnow
%     pause(0.01)
end
%%
x1=13;
x2=112;
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
    w=squeeze(A(n,120:230,:));
    w=mean(w,1);
    a=w<1290;
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
ax_min=TR(:,1)>w+1;
ax_max=TR(:,1)<max(TR(:,1))-w-2;
% ax_max=TR(:,1)<36;

ay_min=TR(:,2)>w+2;
ay_max=TR(:,2)<max(TR(:,2))-w-2;
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
        at=tt>3; % difference between y coordinate
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
yi=[135,143,154,170];
xi=[487.7,546.5,611.6,707];
x0 = 60;
fy = [650:10:750];%for xtick in plotting
fy2= [650:1:750];%for interpolation
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

