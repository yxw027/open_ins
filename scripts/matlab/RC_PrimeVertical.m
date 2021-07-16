% function to compute the radius of curvature in the prime vertical
% by Eun-Hwan Shin, March 2001
%--------------
% Usage: N = rc_PrimeVertical(a, b, lat)
%    - a = semi-major axis of a ellipsoid
%    - e2 = linear eccentricity squared
%    - lat = latitude
function N = rc_PrimeVertical(a, e2, lat);
N = a/sqrt(1-e2 * sin(lat)^2);
