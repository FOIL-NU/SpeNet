l = size (spe,1);
cfit = zeros (l,1);
fy2 = [650:2:775];
for n = 1:l
    try
    ft=fit (fy2',spe(n,:)','gauss1');
    cfit (n) = ft.b1;
    catch
%     plot (fy2',spe (n,:));
%     title (n)
%     hold on
%     plot (ft);
%     drawnow
%     pause (0.1)
%     hold off
    end
end
%%
clf

a1 = cfit>600;
a2 = cfit<800;
a = logical (a1.*a2);
cfit2 = cfit (a);
s = s(a);
c = c(a);
% spe = spe (a,:);
no=no(a);
scatter (cfit2,s,'r');
%%
axis ([650 750 100 1000]);
hold on
scatter (c,s,'g');
std_c = std(c);
std_f = std(cfit2);
disp (sprintf('%%std_centroid = %f',std_c));
disp (sprintf('%%std_fit = %f',std_f));
