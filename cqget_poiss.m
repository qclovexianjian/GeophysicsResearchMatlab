function [ poisson ] = cqget_poiss( vp,vs )

poisson = (vp.^2 - 2*vs.^2)./(2*vp.^2 - vs.^2);

end

