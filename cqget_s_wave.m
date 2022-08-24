function [ vs ] = cqget_s_wave( vp, poisson )

vs = sqrt(vp.^2.*(1-2*poisson)./2./(1-poisson));


end

