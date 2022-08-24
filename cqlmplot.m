function cqlmplot( ...
    depthValueCell, ...
    indexCell, ...
    colorCell, ...
    fillSideCell , ...
    hAx, ...
    lowHighBound, ...
    curveParam)
% Plotting for color coded lithology log
% 
% input
% -----
% depthValueCell = {depth, value} for the lithology log
% indexCell      = {indexVector1, indexVector2, ..., indexVectorN} index vectors for
%                  different lithologies
% colorCell      = {color1, color2, ..., colorN} color coding for lithologies
% fillSideList   = {   1  ,    0   ,...,    1  } 1 = fill value to high bound and 0 = fill
%                  value from low bound
% hAx            = axis you want to draw on 
% ******************* default = gca *********************************
% lowHighBound   = low and high bounds for plotting
% ******************* default = [min(value), max(value)] ************
% curveParam     = plotting parameters for curves
% ******************* default = {'color','k','linew',2} *************

if(~exist('hAx','var') || isempty(hAx))
    hAx = gca;
end

depth = depthValueCell{1}; depth = depth(:);
value = depthValueCell{2}; value = value(:);

if ~exist('lowHighBound','var') || isempty(lowHighBound)
    lowHighBound = [min(value),max(value)];
end
if ~exist('curveParam','var') || isempty(curveParam)
    curveParam = {'color','k','linew',2};
end

% plot curve
plot(depth, value, curveParam{:});
ns = length(depth); % total sample number
% loop into lithology
for n = 1:length(indexCell)
    indexCurrent = sort(indexCell{n});
    nv = 1; % vertex number
    kf = 1; % face number
    nbreak = find(diff(indexCurrent)~=1); % find breaks
    grp = zeros(length(nbreak)+1,2);
    grp(:,1) = [indexCurrent(1);indexCurrent(nbreak+1)]; % group startings
    grp(:,2) = [indexCurrent(nbreak);indexCurrent(end)]; % group endings
    
    % initialize face and vertex
    face = nan*ones(size(grp,1),ns);
    vertex = nan*ones(ns,2);
    
    for k = 1:size(grp,1)
        % first, fill value vertexes
        nPoints = grp(k,2) - grp(k,1) + 1;
        vertex(nv:nv+nPoints-1,1) = value(grp(k,1):grp(k,2)); % fill x
        vertex(nv:nv+nPoints-1,2) = depth(grp(k,1):grp(k,2)); % fill y
        % second, add two vertex according to fill side
        if 
    end
end









end

