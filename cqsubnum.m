function cqsubnum( m , n, fontsize,normpos,h)
% Numbering subplot
% m = row number of total subplots
% n = column number of total subplots
% h = handle of the figure
% fontsize = fontsize
% normpos = [x,y] normalized position default is [0.95,1.2] which put it on right-top of
%           each subplot

if ~exist('h','var') || isempty(h)
    h = gcf;
end
if ~exist('fontsize','var') || isempty(fontsize)
    fontsize = 20;
end
figure(h);
for k = 1:m*n
    subplot(m,n,k);
    try
        font = text(normpos(1), normpos(2), ['(',char(k+96),')'], 'units', 'normalized');
    catch
        font = text(0.95, 1.2, ['(',char(k+96),')'], 'units', 'normalized');
    end
    font.FontSize = fontsize;
end



end

