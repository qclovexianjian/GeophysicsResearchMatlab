function [ xt_filter] = cqfktpfilter_polywin( xt, dt, dx,us, f0, f1, tswp, x0,phi_f,gain  )
% Use polygon windowed LMO method to do the de-doppler operation
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
% gain = put gain in the picking window
%
% output
% ------
% xt_filter = xt data after filtering
% fk = fk matrix
% f = f vector cell array
% k = k vector cell array

xvec = ((1:size(xt,2)) - x0)*dx;
t = 0:dt:(size(xt,1)-1)*dt;

if ~exist('gain','var')||isempty(gain)
    gain = 1;
end
figure;
imagesc2(xvec,t,-xt); colormap gray;
cqwva(xt,dt,xvec,[],gain,[],[],[],'hold');
% let user draw polygons
draw_poly = input('Draw Polygons?(y/n)','s');

% the whole data should be firstly processed then fill in the polygon area
 [ xt_filter ] = ...
    cqfkpfilter( xt, dt, dx,us, f0, f1, tswp  );

while strcmp(draw_poly,'y')
    h = impoly;
    BW = createMask(h); % mask binary matrix
    
    [ p ] = cqpickslope( gca,false );
    [ xt_filter_poly ] = ...
    cqfkpfilter( xt.*BW, dt, dx,us, f0, f1, tswp,phi_f,[mean(p),x0]  );
    
    xt_filter(BW==1) = xt_filter_poly(BW==1);

    draw_poly = input('Draw Polygons?(y/n)','s');
end


end


