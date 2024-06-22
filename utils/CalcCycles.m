function [dila, cons] = CalcCycles(z_bdpupil)
% Input: phase angle of pupil size
% Output: a vector of dilation or constriction events (1 and -1)
z_bdpupil_notnan = z_bdpupil(find(~isnan(z_bdpupil)));
cyc = zeros(size(z_bdpupil_notnan));

dila_angle = angle(hilbert(z_bdpupil_notnan));
diff_dila_angle = [diff(dila_angle); 0];
cyc(find(diff_dila_angle<=-5)) = 1;

cons_angle = angle(hilbert(-z_bdpupil_notnan));
diff_cons_angle = [diff(cons_angle); 0];
cyc(find(diff_cons_angle<=-5)) = -1;
if ~isnan(z_bdpupil(1))
    dila = find(cyc==1);
    cons = find(cyc==-1);
else
    notnan1 = find(~isnan(diff(z_bdpupil)));
    notnan1 = notnan1(1);
    nanlast = notnan1-1;
    dila = find(cyc==1)+nanlast;
    cons = find(cyc==-1)+nanlast;
end
return