function [x50,thickness] = front_metrics(x,profile)
%FRONT_METRICS Center and positive 10--90 width of a left-high front.
z = (profile-min(profile))/(max(profile)-min(profile));
if any(~isfinite(z)); x50=NaN; thickness=NaN; return; end
x90 = crossing(0.9); x50 = crossing(0.5); x10 = crossing(0.1);
thickness = x10-x90;
    function xc = crossing(level)
        j = find(z<=level,1,'first');
        if isempty(j) || j==1; xc=NaN; return; end
        xc = interp1(z(j-1:j),x(j-1:j),level);
    end
end
