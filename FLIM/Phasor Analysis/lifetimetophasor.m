function [g, s] = lifetimetophasor(tau, w)
phi = atan(w*tau);
m = 1./sqrt(1+(w*tau).^2);
[g, s] = pol2cart(phi, m);