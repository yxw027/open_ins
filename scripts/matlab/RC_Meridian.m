% Function to compute the radius of curvature in the meridian
%   by Eun-Hwan Shin, April 2001
%---------------
% Usage: M = rc_meridian(a, e2, lat)
%    - a = semi-major axis
%    - e2 = linear eccentricity squared
%    - lat = latitude
%
function M = rc_meridian(a, e2, lat)
M = a*(1-e2)/(1-e2 * sin(lat)^2)^(3/2);