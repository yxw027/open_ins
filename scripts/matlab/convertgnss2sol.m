%Convert GNSS fromat to be compare

function data_sol = convertgnss2sol(gnss)
deg2rad = 0.017453292519943;
num =length(gnss) ;
column = size(gnss,2);
if(column < 21)
gnss(:,column+1:21) = zeros(num,21 -column);
end
%gnss format
%/* post gnss format*/
%/*GNSS week £¬GNSS second£¨s£©£¬ reserve£¬latitude£¨deg£©£¬longitude£¨deg£©£¬altitude£¬latitude_std£¬longitude_std£¬altitude_std
%north_velocity£¬east_velocity£¬up_velocity£¬north_velocity_std£¬east_velocity_st£¬up_velocity_std, heading,heading_std,numSatellites£¬HDOP£¬sol_age£¬gpsFixType*/
data_sol(:,1) =  gnss(:,2);
data_sol(:,2) =  gnss(:,3) * deg2rad;
data_sol(:,3) =  gnss(:,4) * deg2rad;
data_sol(:,4) =  gnss(:,5);
data_sol(:,5) =  ones(num,1) * 0 ;
data_sol(:,6) =  ones(num,1) * 0 ;
data_sol(:,7) =  ones(num,1) * 0 ;
data_sol(:,8) =  ones(num,1) * 0 ;
data_sol(:,9) =  ones(num,1) * 0 ;
data_sol(:,10) =  ones(num,1) * 0 ;
data_sol(:,11) = gnss(:,9);
end