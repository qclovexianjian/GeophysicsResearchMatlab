function [ xt_filter,fk, f, k,pwin  ] = cqfktpfilter_window( xt, dt, dx,us, f0, f1, tswp, x0,pshift_window,margin,gain  )
% Use windowed LMO method to do the de-doppler operation
% This function is a higher-level function of cqfktpfilter
%
% input
% -----
% xt = xt domain data, either cdp, crg or csg
% dt = time sampling rate
% dx = space sampling rate
% us = boat speed
% f0 = start sweeping frequency
% f1 = end sweeping frequency (assuming linear sweeping)
% tswp = sweeping duration in seconds
% x0 = index where offset = 0
% pshift_window = [pshift_1, window_1_xleft, window_1_xright, window_1_ytop, window_1_ybottom]
%                 [pshift_2, window_2_xleft, window_2_xright, window_2_ytop, window_2_ybottom]
%                               ...
%                 [pshift_N, window_N_xleft, window_N_xright, window_N_ytop, window_N_ybottom]
% margin = [x_margin, y_margin] indicades margin area in x and y direction
% gain = put gain in the picking window
%
% output
% ------
% xt_filter = xt data after filtering
% fk = fk matrix 
% f = f vector cell array
% k = k vector cell array
% pwin = dividing window parameters stored for further usage

% if user dosen't provide pshift_window, then we will use an interactive program to let
% user devide the data into regions
xvec = ((1:size(xt,2)) - x0)*dx;
t = 0:dt:(size(xt,1)-1)*dt;
if ~exist('pshift_window','var')||isempty(pshift_window)
    if ~exist('gain','var')||isempty(gain)
        gain = 1;
    end
    figure;
    cqwva(xt,dt,xvec,[],gain); hold on;
    % let user draw horizontal divide
    title('Draw Horizontal Divide')
    add = input('add more horizontal divide?(y/n)','s');
    h = {}; n = 1; ydiv = []; xdiv = [];
    while strcmp(add,'y')
        h{n} = ginput2(1);
        plot([xvec(1),xvec(end)],[h{n}(2),h{n}(2)],'r','linew',3);
        add = input('add more horizontal divide?(y/n)','s');
        n = n + 1;
    end
    % let user draw vertical divide
    title('Draw Vertical Divide')
    add = input('add more vertical divide?(y/n)','s');
    ver = {}; n = 1;
    while strcmp(add,'y')
        ver{n} = ginput2(1);
        plot([ver{n}(1),ver{n}(1)],[t(1),t(end)],'b','linew',3);
        add = input('add more vertical divide?(y/n)','s');
        n = n + 1;
    end
    % let user choose slope in each region
    nwin = (length(h)+1) * (length(ver)+1);
    pshift_window = zeros(nwin,5);
    for m=1:length(ver) 
        xdiv(m) = ver{m}(1); 
    end
    xdiv = [xvec(1),xdiv,xvec(end)];
    for m=1:length(h) 
        ydiv(m) = h{m}(2); 
    end
    ydiv = [t(1),ydiv,t(end)];
    for m = 1:(nwin)
        nxdiv = ceil(m/(length(h)+1));
        nydiv = mod(m,(length(h)+1));
        if nydiv == 0
            nydiv = length(h)+1;
        end
        [ this_p ] = cqpickslope( gca );
        if isempty(this_p)
            this_p = 0;
        else
            this_p = this_p(end);
        end
        pshift_window(m,:) = [this_p,xdiv(nxdiv),xdiv(nxdiv+1),ydiv(nydiv),ydiv(nydiv+1)];
    end
end

xt_filter = zeros(size(xt));
pwin = pshift_window; % simplify the notation

fk = cell(size(pwin,1),1);
f = fk;
k = fk;

for n = 1:size(pwin,1)
    nx = close_pos(xvec,pwin(n,[2,3]));
    ny = close_pos(t,pwin(n,[4,5]));
    nx_pad = close_pos(xvec,pwin(n,[2,3])+[-margin(1),margin(1)]);
    ny_pad = close_pos(t,pwin(n,[4,5])+[-margin(2),margin(2)]);
    relative_x = nx - nx_pad(1) + 1; % relative position of x in nx_pad
    relative_y = ny - ny_pad(1) + 1; % relative position of y in ny_pad
    xt_temp = xt(ny_pad(1):ny_pad(2),nx_pad(1):nx_pad(2));
    [ xt_filter_temp,~,fk{n},f{n},k{n} ] = ...
        cqfkpfilter( xt_temp, dt, dx,us, f0, f1, tswp, [pwin(n,1),x0]  );
    xt_filter_temp = xt_filter_temp(relative_y(1):relative_y(2),relative_x(1):relative_x(2));
    xt_filter(ny(1):ny(2),nx(1):nx(2)) = xt_filter_temp;
end

end

function npos = close_pos(xpos,pos)
% return the nearest index of the position in xpos
% mg is margin
npos = zeros(size(pos));
for m = 1:length(pos)
    [~,npos(m)] = min(abs(xpos-pos(m)));
end
end

