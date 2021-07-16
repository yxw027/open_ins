% Convert Euler angles to Dirction Cosine Matrix
% function dc = euler2dcm(roll, pitch, heading)
% Input: three Euler angles in rad.

function dc = euler2dcm(roll, pitch, heading)

dc = zeros(3);
cr = cos(roll); cp = cos(pitch); ch = cos(heading);
sr = sin(roll); sp = sin(pitch); sh = sin(heading);

dc(1,1) = cp * ch ;
dc(1,2) = -cr*sh + sr*sp*ch;
dc(1,3) = sr*sh + cr*sp*ch ;

dc(2,1) = cp * sh;
dc(2,2) = cr*ch + sr*sp*sh;
dc(2,3) = -sr * ch + cr * sp * sh;

dc(3,1) = - sp;
dc(3,2) = sr * cp;
dc(3,3) = cr * cp;