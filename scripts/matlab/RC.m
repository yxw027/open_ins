% Function to compute the radius of curvature in the meridian and the prime
% vertical
% function [M,N] = RC(a, e2, lat)
%    - a = semi-major axis
%    - e2 = linear eccentricity squared
%    - lat = latitude

function [M,N] = RC(a, e2, lat)
M = a*(1-e2)/(1-e2 * sin(lat)^2)^(3/2);
N = a/sqrt(1-e2 * sin(lat)^2);