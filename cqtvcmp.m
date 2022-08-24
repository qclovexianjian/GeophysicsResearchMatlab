function [ d_comp ] = cqtvcmp( d, tori, tdelay )
% Time varying compress d


tori = tori(:);
tdelay = tdelay(:);
d_comp = [];
tq = tori+tdelay;
for k = 1:size(d,2)
    d_comp(:,k) = interp1(tori,d(:,k),tq,'linear',d(end));
end


end

