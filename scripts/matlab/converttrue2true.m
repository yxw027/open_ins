function data_true = converttrue2true(true,offset)
deg2rad = 0.017453292519943;
data_true = true(:,2:end);
data_true(:,2:3) = data_true(:,2:3) * deg2rad;
data_true(:,8:10) = data_true(:,8:10) * deg2rad;


WGS84=GetWGS84;
M=RC_Meridian(WGS84.a, WGS84.e2, data_true(1,2));
N=RC_PrimeVertical(WGS84.a, WGS84.e2, data_true(1,2));
for i=1:length(data_true)
        C_bn = euler2dcm(data_true(i,8), data_true(i,9), data_true(i,10));
        la_r = C_bn * offset;
        d_l(1) = la_r(1) / (M);
        d_l(2) = la_r(2) / (N * cos(data_true(i,2))) ;
        d_l(3) = -la_r(3);
        for j =2:4
        data_true(i,j) = data_true(i,j) +d_l(j-1);
        end    
 end

end