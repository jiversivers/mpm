function [a1, a2, a2Percent] = getPhasorFractions(gp, sp, tau1, tau2, w)

[g1, s1] = lifetimetophasor(tau1, w);
[g2, s2] = lifetimetophasor(tau2, w);

a1 = sqrt((gp-g2).^2 + (sp-s2).^2);
a2 = sqrt((gp-g1).^2 + (sp-s1).^2);
a2Percent = a2./(sqrt((g2-g1)^2+(s2-s1)^2));