%Convert UBLOX fromat to be compare

function data_sol = convertublox2sol(ublox)
deg2rad = 0.017453292519943;
num =length(ublox) ;
column = sizeof(ublox,2);
if(column < 24)
ublox(:,column+1:24) = zeros(num,24 -column);
end
%ublox format
%/* ublox csv format*/
data_sol(:,1) = ublox(:,2)/1000;
data_sol(:,2:3) = ublox(:,9:10) * deg2rad;
data_sol(:,4) = ublox(:,11);
data_sol(:,5) = ublox(:,19);
data_sol(:,6) = ublox(:,18);
data_sol(:,7) = ublox(:,20);
data_sol(:,8) = ones(size(ublox,1),1) * 0;
data_sol(:,9) =nes(size(ublox,1),1) * 0;
data_sol(:,10) = ublox(:,24) * deg2rad;
end