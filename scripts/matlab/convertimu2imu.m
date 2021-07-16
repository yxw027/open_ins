function data_imu = convertimu2imu(imu)
%data format 
%Gnss_week  Gnss_second reset  Acc_X  Acc_Y  Acc_Z  Gyro_X  Gyro_Y  Gyro_Z
%    ()       (s)        (m/s^2)(m/s^2)(m/s^2)   (deg/s) (deg/s) (deg/s)
%IMU format
%data format 
%Gnss_second   Gyro_X  Gyro_Y  Gyro_Z Acc_X  Acc_Y  Acc_Z 
%    (s)     (rad/s) (rad/s) (rad/s)(m/s^2)(m/s^2)(m/s^2)
deg2rad = 0.017453292519943;
%num =length(sol) ;
%column = sizeof(fil,2);
% if(column < 21)
% fil(:,column+1:21) = zeros(num,21 -column);
% end
data_imu(:,1) = imu(:,2); 
data_imu(:,2:4) = imu(:,7:9) *deg2rad;
data_imu(:,5:7) = imu(:,4:6);
end