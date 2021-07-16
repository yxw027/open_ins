function   ret = show_result_all(n_sol, imufilename,gnssfilename,solfilename,outputfilefolder,filfilename,reffilename,cfg_info)
ret=1 ;

mkdir(outputfilefolder);
debug_filename = [outputfilefolder 'debug.txt'];
fp_debug = fopen(debug_filename,'wt');

if ~isempty(imufilename)
    rawimu = load(imufilename);
    imu = convertimu2imu(rawimu);
    clear rawimu;
    
    mkdir([outputfilefolder 'imu/']);
    Show_imu(imu,[outputfilefolder 'imu/'], fp_debug);
end

rawgnss = load(gnssfilename);
rawsol = load(solfilename);
%rawfil = load(filfilename);
rawref = load(reffilename);
%数据转化

ref = converttrue2true(rawref,(cfg_info.leverarm.ref2gnss));
clear rawref;
gnss_sol= convertgnss2sol(rawgnss);
clear rawgnss;
sol = convertsol2sol(rawsol);
clear rawsol;
% fil = convertfil2fil(rawfil);

%需要完成分析的工作
%1、绘制IMU的曲线，输出IMU起始时间，终止时间，采样频率，是否地丢数，频率是否稳定

%分析error type= 1,分析GNSS误差 保存误差，统计结果，保存误差曲线
mkdir([outputfilefolder 'gnss/']);
% Show_sol(1,gnss_sol,ref,cfg_info.dsc_nogs,cfg_info.t_nogps,[outputfilefolder 'gnss/']);

mkdir([outputfilefolder 'ins/']);
Show_sol(4, n_sol, sol,ref,cfg_info.dsc_nogs,cfg_info.t_nogps, ...
    cfg_info.cdf_fh,cfg_info.cep_fh, cfg_info.legend,...
    [outputfilefolder 'ins/']);

fclose(fp_debug);
end
