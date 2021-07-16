%
% function to get WGS84 coefficients
% by Eun-Hwan Shin, 2001
%  function wgs84 = GetWGS84
%
function wgs84 = GetWGS84
wgs84.a = 6378137.0;
wgs84.b = 6356752.3142;
wgs84.f = 1.0/298.2572236;
wgs84.we = 7.292115147e-5;
wgs84.e2 = 1 - (wgs84.b/wgs84.a)^2;
wgs84.gm = 3986004.418e+8;
wgs84.j2 = 1.082626683e-3;
wgs84.j3 = -2.5327e-6;
