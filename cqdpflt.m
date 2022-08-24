function [ gather ] = cqdpflt( xs0, xr0, us, ur, z,xref,theta, ...
    v, f0, f1, tswp, taper, tuncor, tcor, dt,use_chirp, v_rho,custom_pilot )
% Model Doppler effects for a flat reflector as well as a dipping reflector
% 
% input
% -----
% xs0 = vector for source locations
% xr0 = vector for receiver locations
% *****
% program determine gather type based on xs0 and xr0 size
% If ns = 1, program assumes to generate common-shot-gather
% If nr = 1, program assumes to generate common-receiver-gather
% If ns = nr, program assumes to generate common-midpoint-gather
% *****
% us = scalor for source speed (m/s)
% ur = scalor for receiver speed. Both speeds are constants. (m/s)
% z = vector for reflector depth 
% xref = reference x for the depth. it makes no difference for flat reflector. but is
%        required for dipping reflection.
% theta = vector of size(z) to indicate degree angles for the reflectors. Note, positive
%         theta is down-dipping and negative theta is up-dipping
% v = scalor for medium (usually) water speed
% f0 = starting sweeping frequency in Hz
% f1 = ending sweeping frequency in Hz
% tswp = sweeping time length in seconds
% taper = taper at both ends in seconds
% tuncor = maximum time in uncorrelated gather in seconds
% tcor = maximum time in correlated gather in seconds
% dt = sampling rate in seconds
% use_chirp = if empty or not given, program use cqchirp
%             else user should give one of the followings:
%               'linear','quadratic', and 'logarithmic' OR
%               'custom'
%              and program will use matlab default chirp function
% v_rho = (optional) AVO effects would be added into solution when it is supplied.
%         v_rho should be [vp1,vs1,rho1,vp2,vs2,rho2]
% custom_pilot = if user give 'custom' for use_chirp parameter, he has to give the
%         custom_pilot
%
% output
% ------
% gather = structure with output information
%    dcor = correlated gather
%    duncor = uncorrelated gather


if ~exist('use_chirp','var') || isempty(use_chirp)
    if_use_chirp = false;
else
    if_use_chirp = true;
    if strcmp(use_chirp,'custom')
        t_pilot = 0:dt:(length(custom_pilot)-1)*dt;
        t_pilot = t_pilot(:); % this is for custom pilot trace
    end
end
% check input
theta = theta/180*pi;
ns = length(xs0);
nr = length(xr0);
if ns~=nr && nr~=1 && ns~=1
    error('Invalid xs0 and xr0 length!');
end

t = 0:dt:tuncor; % get recording time 
t = t(:);
nsamp = length(t); % total samples per trace
tuncor = t(end); % in case tuncor is not an integer times dt

% Create pilot trace
t_pilot = 0:dt:tswp; t_pilot = t_pilot(:);

if ~if_use_chirp
    tr_plt = cqchirp(t_pilot,f0,tswp,f1,0,taper);
else
    if strcmp(use_chirp,'custom')
        tr_plt = custom_pilot(:);
    else
        tr_plt = chirp(t_pilot,f0,tswp,f1,use_chirp);
    end
end

t_align = cell(length(z),1);

if ns == nr
    dcor = zeros(nsamp,ns);
    gather.type = 'common-midpoint-gather';
elseif ns == 1
    dcor = zeros(nsamp,nr);
    gather.type = 'common-shot-gather';
else
    dcor = zeros(nsamp,ns);
    gather.type = 'common-receiver-gather';
end
duncor = dcor;
gather.inc_ang = {};
gather.coef = {};



for nz = 1:length(z)
if ns == nr
    % create common-midpoint-gather/random-shoot
    for k = 1:ns
        [ts,gather.inc_ang{nz,k},gather.coef{nz,k}] = search(xs0(k),xr0(k),z(nz),nz);
        % create uncorrelated records
        if ~if_use_chirp
            duncor(:,k) = duncor(:,k) + gather.coef{nz,k}.*cqchirp(ts,f0,tswp,f1,0,taper);
        else
            try
                duncor(:,k) = duncor(:,k) + gather.coef{nz,k}.*chirp(ts,f0,tswp,f1,use_chirp);
            catch
                if strcmp(use_chirp,'custom')
                    duncor(:,k) = duncor(:,k) + gather.coef{nz,k}.*...
                        interp1(t_pilot,custom_pilot,ts,'spline',0);
                else
                    error('Not a valid chirp name! Correct use_chirp argument!');
                end
            end
        end
        t_align{nz}(k) = interp1(ts, t, 0, 'linear', 'extrap');
    end
    
elseif ns == 1
    % create common-shot-gather
    for k = 1:nr
        [ts,gather.inc_ang{nz,k},gather.coef{nz,k}] = search(xs0,xr0(k),z(nz),nz);
        % create uncorrelated records
        if ~if_use_chirp
            duncor(:,k) = duncor(:,k) + gather.coef{nz,k}.*cqchirp(ts,f0,tswp,f1,0,taper);
        else
            try
                duncor(:,k) = duncor(:,k) + gather.coef{nz,k}.*chirp(ts,f0,tswp,f1,use_chirp);
            catch
                if strcmp(use_chirp,'custom')
                    duncor(:,k) = duncor(:,k) + gather.coef{nz,k}.*...
                        interp1(t_pilot,custom_pilot,ts,'spline',0);
                else
                    error('Not a valid chirp name! Correct use_chirp argument!');
                end
            end
        end
        t_align{nz}(k) = interp1(ts, t, 0, 'linear', 'extrap');
    end
    
else
    for k = 1:ns
        [ts,gather.inc_ang{nz,k},gather.coef{nz,k}] = search(xs0(k),xr0,z(nz),nz);
        % create uncorrelated records
        if ~if_use_chirp
            duncor(:,k) = duncor(:,k) + gather.coef{nz,k}.*cqchirp(ts,f0,tswp,f1,0,taper);
        else
            try
                duncor(:,k) = duncor(:,k) + gather.coef{nz,k}.*chirp(ts,f0,tswp,f1,use_chirp);
            catch
                if strcmp(use_chirp,'custom')
                    duncor(:,k) = duncor(:,k) + gather.coef{nz,k}.*...
                        interp1(t_pilot,custom_pilot,ts,'spline',0);
                else
                    error('Not a valid chirp name! Correct use_chirp argument!');
                end
            end
        end
        t_align{nz}(k) = interp1(ts, t, 0, 'linear', 'extrap');
    end
end
end

for ncorr = 1:size(duncor,2)
    % create correlated records
    [s,lag] = xcorr(duncor(:,ncorr),tr_plt);
    nt0 = find(lag==0);
    dcor(:,ncorr) = s(nt0:nt0+nsamp-1);
end

dcor = dcor(1:floor(tcor/dt),:);


gather.data_correlate = dcor;
gather.data_uncorrelate = duncor;
gather.pilot_trace = tr_plt;
gather.dt = dt;
gather.source_speed = us;
gather.receiver_speed = ur;
gather.source_x = xs0;
gather.receiver_x = xr0;
gather.freq_low = f0;
gather.freq_high = f1;
gather.taper = taper;
gather.sweeping_time = tswp;
gather.reflector_depth = z;
gather.nmo_time = t_align;
gather.tuncor = tuncor;
gather.tcor = tcor;



% =================
    function [ts,inc_ang,coef] = search(xs0, xr0,z,nrlt)
        % search for source time
        % xs0 and xr0 are only scalar
        % nrlt represents which reflector you are working on
        if theta(nrlt)==0
            D = - xr0 - t*ur + xs0 + t*us; % D is broadcast by tuncor
            A = us^2 - v^2;
            B = -2 * D * us; % broadcast
            C = D.^2 + 4*z^2; % broadcast
        else
            x_theta = xref(nrlt) - z*cot(theta(nrlt));
            D = xs0 + t*us - xr0 - t*ur;
            H = xs0 + us*t - x_theta;
            A = us^2 - v^2;
            B = -8 * H * us * sin(theta(nrlt))^2 - 2 * D * us + ...
                us * 4 * sin(theta(nrlt))^2*(D + H);
            C = 4 * sin(theta(nrlt))^2 * ( H.^2 - D.*H ) + D.^2;
        end
        delt = (-B - sqrt(B.^2-4*A*C))/(2*A); % time from source to receiver
        ts = t - delt;
        % calculate incident angle
        xs = xs0 + us * ts;
        if theta(nrlt)==0
            z0 = z;
        else
            z0 = (xs - x_theta) * sin(theta(nrlt)); % broad cast z0
        end        
        E = xs0 - xr0 + us*delt;
        F = delt * v;
        inc_ang = acos((4*z0.^2+F.^2-E.^2)/4./F./z0);
        if exist('v_rho','var')&&~isempty(v_rho)
            coef=zoeppritz(v_rho(3),v_rho(1),v_rho(2),v_rho(6),...
                v_rho(4),v_rho(5),1,1,0,inc_ang/pi*180); % coef = real + i*imag
            coef(imag(coef)~=0) = 0; % ignore angles after critical angle
        else
            coef = ones(size(inc_ang));
        end
        coef = coef(:);
    end


end

