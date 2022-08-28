% Load spectrum file
clear
clc;

fnamein=sprintf('cali.nd2');
A = SPLMload(fnamein,'nd2','double');
A = rot90(A,3);
% A = rot90(A,3);
% A = rot90(A,3);
figure(1);imagesc(squeeze(mean(A,3))); axis off; axis image; axis image; colormap(gray)
A=double(squeeze(A));
A= A - 200;

cali=mean(A,3);
cali=mean(cali,1);
plot(cali)

%Calibration data
yi=[219,228,241];
xi=[487.7,546.5,611.6];
res=(xi(end)-xi(1))/(yi(end)-yi(1))
