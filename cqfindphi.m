function [ phi,f ] = cqfindphi( xt, xt0, dt, p, us, f0, f1 )
% The function finds normalized phi for doppler factor = 1
% 
% input
% -------
% xt = one trace with phase distortion
% xt0 = one trace at the same location without phase distortion
% Note: both xt and xt0 should be uncorrelated and should be only containing first
% arrival. xt and xt0 should be in the same length
% dt = sampling rate in seconds
% p = slope
% us = boat speed
% f0 = start frequency
% f1 = end freuqnecy 
%
% output
% ------
% phi = phase distortion function by Doppler Effects
% f = corresponding frequnecy vector


% FFT xt and xt0 to frequency domain and set 0 for frequency out of f0 ~ f1
[XT, f] = cqfft(xt,dt,0);
[XT0] = cqfft(xt0,dt,0);
XT(f<f0 | f>f1) = 0;
XT0(f<f0 | f>f1) = 0;

% Compute the subtraction of their phase spectra
phi = unwrap(unwrap(angle(XT)) - unwrap(angle(XT0)));

% Normalize the phase spectra by its doppler factor which should be calculated as -p*us
phi = phi./(-p*us);
phi(f<f0 | f>f1) = 0;


end

