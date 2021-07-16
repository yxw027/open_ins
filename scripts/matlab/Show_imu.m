% Program to read the IMU data file, adjust the axis, display and save it.
% The output is the median of the sampling interval.
% function result=show_IMU(fname)
% By Xiaoji Niu, Mar 2007

function mean_dt=Show_imu(xx,outputfilefolder,fp_debug)
%fname='dr_static1.imu';
% [N,~]=size(xx);

idx = find(xx(:,1));
N = length(idx); 
t = xx(idx, 1);
xx = xx(idx, :);

dt=diff(t); mean_dt=median(dt);
% Brief description
fprintf(fp_debug,'IMU start time:%10.4f (second),end time:%10.4f (second),interval:%5.3f (second) \n',xx(1,1),xx(end,1),mean_dt);
% Exception description
for i = 1:N-1
    if (dt(i) > 2 * mean_dt || dt(i) < 0.5 * mean_dt)
        fprintf(fp_debug,'IMU Abnormal time interval: time:%10.4f (second),  interval:%10.4f (second)\n',xx(i+1,1),dt(i));
    end
end

% Plot the signals to check
t0=floor(xx(1,1)/10000)*10000;
fh = figure('Visible','on');
ax(1)=subplot(311); plot(xx(1:N,1)-t0,xx(1:N,2:4)*180/pi), ylabel('Gyro (deg/s)'), %./(dt*ones(1,3))
grid on; legend('Forward','Right','Down'),
if t0==0
    xlabel(' Time (sec)');
else
    xlabel(['Time -' int2str(t0) ' (sec)']);
end
ax(2)=subplot(312); plot(xx(1:N,1)-t0,xx(1:N,5:7)), ylabel('Acc (m/s^2)'), %./(dt*ones(1,3))
grid on; legend('Forward','Right','Down'),
if t0==0
    xlabel(' Time (sec)');
else
    xlabel(['Time -' int2str(t0) ' (sec)']);
end
ax(3)=subplot(313); plot(xx(2:N,1)-t0,dt), ylabel('Interval (sec)'), grid on;
if t0==0
    xlabel(' Time (sec)');
else
    xlabel(['Time -' int2str(t0) ' (sec)']);
end
linkaxes(ax,'x');
title("IMU");
saveas(fh,[outputfilefolder 'imu.fig']);
dt = sort(abs(dt-mean_dt)/mean_dt);
err_95p =dt(round(0.95*(N-1)));
err_99p =dt(round(0.99*(N-1)));
err_xp =dt(round(0.9997*(N-1)));
err_xxp =dt(round(0.999997*(N-1)));
fprintf(fp_debug,'IMU Abnormal time stability :%10.4f,%10.4f,%10.4f,%10.4f \n',err_95p,err_99p,err_xp,err_xxp);

% figure, % plot the histogram of the sampling interval.
%  hist(dt, 0.005:0.001:ceil(max(dt)/0.1)*0.1); grid on;
%  xlabel('Interval (sec)'), ylabel('Point number'), title('Histograph of the DR sampling interval');
%  title("IMU interval");
% saveas(gcf,[outputfilefolder 'imuinterval.fig']);
end
