function [ filter,amp,f ] = cqgenwav( fq,val,dt,n, interp_method )
% Frequency domain genwav equation similar to Dr. Hilterman's
% Input
% -----
% fq = query frequency in Hz
% val = amplitude for the query frequency fq (val length should equal to fq length)
% dt = sampling rate in second
% n = sample number in the total filter (need to be an odd number)
% interp_method = valid interpolation method allowed, refer to interp built-in function
%
% Output
% ------
% filter = the designed filter with time zero at middle
% amp = designed amplitude spectrum
% f = frequency axis for amp

df = 1/dt/n;
f = 0:df:1/2/dt;
f = f(:);
fq = fq(:);
val = val(:);

amp = interp1(fq,val,f,interp_method,0);
amp_ifft = [amp;flipud(amp(2:end-1))];

filter = [ifftshift(real(ifft(amp_ifft)));0];


end

