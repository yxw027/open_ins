function data_gnss = convertgnss2gnss(gnss)
deg2rad = 0.017453292519943;
%num =length(sol) ;
%column = sizeof(fil,2);
% if(column < 21)
% fil(:,column+1:21) = zeros(num,21 -column);
% end
data_gnss = gnss(:,2:end);
data_gnss(:,2:3) = gnss(:,2:3) *deg2rad;
end