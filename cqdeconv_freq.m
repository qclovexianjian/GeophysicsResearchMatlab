function [ d_decon ] = cqdeconv_freq( d, WL,dt,fWL, prewt )
% The function apply freq domain decon assuming knowing wavelets spectrum
% input
% -----
% d = 2D matrix of seismic data,column vector if trace number == 1
% WL = the vector representing wavelet's spectrum. 
% dt = sampling rate in seconds
% fWL = sampling frequency vector for WL vector
% prewt = prewhitenning (max(abs(WL))*0.01)
% 
% output
% ------
% d_decon = deconvolution result of the input seismic data

if ~exist('prewt','var')||isempty(prewt)
    prewt = 0.1;
end
prewt = prewt * max(abs(WL));
invWL = 1./(WL + prewt); % spectrum of the inverse filter after prewhitening
[D,fD] = cqfft(d,dt,0);
[ invWL ] = cqcomplex_interp( invWL, fWL, fD, 12 ); % Complex number interpolation
invWL = invWL(:);
invWLmat = repmat(invWL,1,size(D,2));
D_DECON = D.*invWLmat;
D_DECON = [real(D_DECON);flipud(real(D_DECON(2:end-1,:)))] + ...
    1i* [imag(D_DECON);-flipud(imag(D_DECON(2:end-1,:)))];
d_decon = real(ifft(D_DECON,[],1));

end

