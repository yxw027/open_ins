% function to calculate mean, rms, and max of an error vector.

function stat=sum_statistic(error)
stat(1)=mean(abs(error));
stat(2)=norm(error)/sqrt(length(error));
stat(3)=max(abs(error));
