%Convert GNSS fromat to be compare

function data_sol = convertsol2sol(sol)
deg2rad = 0.017453292519943;
%num =length(sol) ;
%column = sizeof(sol,2);
% if(column < 21)
% sol(:,column+1:21) = zeros(num,21 -column);
% end
%gnss format
%/* post gnss format*/
%/*GNSS week £¬GNSS second£¨s£©£¬ reserve£¬latitude£¨deg£©£¬longitude£¨deg£©£¬altitude£¬latitude_std£¬longitude_std£¬altitude_std
%north_velocity£¬east_velocity£¬up_velocity£¬north_velocity_std£¬east_velocity_st£¬up_velocity_std, heading,heading_std,numSatellites£¬HDOP£¬sol_age£¬gpsFixType*/
data_sol = sol(:,2:end);
data_sol(:,2:3) = data_sol(:,2:3) *deg2rad;
data_sol(:,8:10) = data_sol(:,8:10) * deg2rad;
data_sol(:,11:13) = data_sol(:,11:13) *deg2rad /3600;
data_sol(:,14:16) = data_sol(:,14:16) /(10e5);
end