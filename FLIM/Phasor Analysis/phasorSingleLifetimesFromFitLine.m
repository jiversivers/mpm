function [tau1, tau2] = phasorSingleLifetimesFromFitLine(varargin)

% f is a linear fit object, m and b are slope and y-int respectively, w is
% the angular fredquency of the laser
if nargin == 2
    f = varargin{1};
    w = varargin{2};
    m = f.p1;
    b = f.p2;
elseif nargin == 3
    m = varargin{1};
    b = varargin{2};
    w = varargin{3};
end

% Equation for line and circle solved for point(s) of equality
x = (-(2*m*b-1)+[-1 1]*sqrt((2*m*b-1)^2-4*(m^2+1)*b^2))/(2*(m^2+1));
y = m*x+b;

tau = phasortolifetime(x, y, w);
tau = sort(tau);
tau1 = tau(1);
tau2 = tau(2);

if ~isreal(tau1) || ~isreal(tau2)
    warning(sprintf("Line does not intersect universal circle. \nThis may happen in cases of single-lifetimes resulting in a 'near-miss.'\nIn this case, the returned approxiamtion may be acceptable for further anaylsis."))
    tau1 = real(tau1);
    tau2 = real(tau2);
end