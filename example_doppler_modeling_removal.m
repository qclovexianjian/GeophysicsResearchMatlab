% Test the dipping reflector modeling program for Doppler effects
xs0 = 600;
xr0 = 300;
xref = 0;
theta = 30;
us = 0;
ur = 0;
z = 0;
v = 1000;
f0 = 18;
f1 = 75;
tswp = 18;
taper = 0.5;
tuncor = 20;
tcor = 2;
dt = 0.001;

% model 1
[ gather ] = cqdpflt( xs0, xr0, us, ur, z,xref,theta, ...
    v, f0, f1, tswp, taper, tuncor, tcor, dt );
% Hand calculation of the result
zs = xs0*sin(theta/180*pi);
d = sqrt((xs0-xr0)^2 + 4*zs^2 - 4*(xs0-xr0)*zs*sin(theta/180*pi));
t_reflt = d/v;
cqwva(gather.data_correlate,dt); title(['Peak should arrive at ',num2str(t_reflt)]);

% model 2
xref = 1200;
theta = -30;
[ gather ] = cqdpflt( xs0, xr0, us, ur, z,xref,theta, ...
    v, f0, f1, tswp, taper, tuncor, tcor, dt );

zs = (xs0-xref)*sin(theta/180*pi);
d = sqrt((xs0-xr0)^2 + 4*zs^2 - 4*(xs0-xr0)*zs*sin(theta/180*pi));
t_reflt = d/v;
cqwva(gather.data_correlate,dt); title(['Peak should arrive at ',num2str(t_reflt)]);

% model 3
xref = 0; theta = 30; xr0 = 600;dx = 5;
xs0 = [600:dx:1200];
[ gather ] = cqdpflt( xs0, xr0, us, ur, z,xref,theta, ...
    v, f0, f1, tswp, taper, tuncor, tcor, dt );
cqwva(gather.data_correlate,dt,xs0); 
title('Down-dip 30 degree at (x = 0, z = 0) receiver at x = 600, boat not moving')

us = -3;
[ gather ] = cqdpflt( xs0, xr0, us, ur, z,xref,theta, ...
    v, f0, f1, tswp, taper, tuncor, tcor, dt );
cqwva(gather.data_correlate,dt,xs0); 
title('Down-dip 30 degree at (x = 0, z = 0) receiver at x = 600, boat moving at 3m/s')

% try to remove that doppler factor
dfilter = cqfkpfilter( gather.data_correlate, dt, dx,us, f0, f1, tswp );
cqwva(dfilter,dt,xs0); 
title('Doppler effect removed for model(with chen method): Down-dip 30 degree at (x = 0, z = 0) receiver at x = 600, boat moving at 3m/s')

duncor = gather.data_uncorrelate;
duncor = cqslantshift(duncor,us,dt,dx);
dcorfilter = cqcorr(duncor,gather.pilot_trace);
dcorfilter = cqslantshift(dcorfilter,-us,dt,dx);
dcorfilter = dcorfilter(1:size(dfilter,1),:);
cqwva(dcorfilter,dt,xs0); 
title('Doppler effect removed for model(with Hampson method): Down-dip 30 degree at (x = 0, z = 0) receiver at x = 600, boat moving at 3m/s')