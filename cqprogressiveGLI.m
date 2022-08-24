function [ inv_imp, f ] = cqprogressiveGLI( d, dt, f1, f2, NfIter, opt )
% The function will run progressive GLI for the input seismogram
%
% input
% -------
% d = input seismic vector
% dt = seismic sampling rate in seconds
% f1, f2 = filter high-end frequency corners. The first low-pass filter will be
% 0-0-f1(1)-f2(2) and the last low-pass filter will be 0-0-f1(2)-f2(2)
% NfIter = number of frequency iterations. The frequency window is
% computed as:
%        [0-0-linspace(f1(1),f1(2),NfIter)-linspace(f2(1),f2(2),NfIter)]
% opt = option structure from glis.opt function
%       The following fields will be abandoned so don't bother with it:
%           1. ShowProgress is set to be false. This program will give you neccessary
%           progress reporting in the command window
%           2. Auto is set to be true for not waiting for user input
%        User is required to give a relatively good source estimation
%        The initial model could be just a constant scalar
%
% output
% -------
% inv_imp = Inverted Impedance Profile in cell arrays
% f = inversion frequency windows in matrix,each row is a frequency window

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% check input and give an option structure for 0 input
if nargin == 0
    inv_imp = glis.opt;
    opt.ShowProgress = false;
    opt.Auto = true;
    return;
end
opt.ShowProgress = false;
opt.Auto = true;
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% calculate frequency inversion window
f = zeros(NfIter,4);
f(:,3) = linspace(f1(1),f1(2),NfIter)';
f(:,4) = linspace(f2(1),f2(2),NfIter)';
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Output some string
disp('+++++++++++++++++++++++++')
disp('Progressive GLI frequency window:')
disp(f)
disp('Start...')
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% some initiation
inv_imp = cell(NfIter,1);
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Iterate for each frequency window
source = opt.Source;
for k = 1:NfIter
    disp(['Inverting window = ',num2str(k)])
    % create filter and convolve with the seismic wavelet as the new wavelet
    [filter] = cqgenwav([0,1,f(k,3),f(k,4)],[1,1,1,0],dt, 2001,'linear');
    opt.Source = cqi_conv(source, filter);
    % filter the input seismic with the low-pass filter
    d_inv = cqi_conv(d, filter);
    % Conventional GLI inversion
    r = glis.dogli(d_inv, opt);
    % Update impedance profile
    inv_imp{k} = r.m;
    opt.InitModel = r.m;
end
disp('+++++++++++++++++++ All Window Finished ++++++++++++++++')



end

