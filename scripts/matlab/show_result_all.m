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
%����ת��

ref = converttrue2true(rawref,(cfg_info.leverarm.ref2gnss));
clear rawref;
gnss_sol= convertgnss2sol(rawgnss);
clear rawgnss;
sol = convertsol2sol(rawsol);
clear rawsol;
% fil = convertfil2fil(rawfil);

%��Ҫ��ɷ����Ĺ���
%1������IMU�����ߣ����IMU��ʼʱ�䣬��ֹʱ�䣬����Ƶ�ʣ��Ƿ�ض�����Ƶ���Ƿ��ȶ�

%����error type= 1,����GNSS��� ������ͳ�ƽ���������������
mkdir([outputfilefolder 'gnss/']);
% Show_sol(1,gnss_sol,ref,cfg_info.dsc_nogs,cfg_info.t_nogps,[outputfilefolder 'gnss/']);

mkdir([outputfilefolder 'ins/']);
Show_sol(4, n_sol, sol,ref,cfg_info.dsc_nogs,cfg_info.t_nogps, ...
    cfg_info.cdf_fh,cfg_info.cep_fh, cfg_info.legend,...
    [outputfilefolder 'ins/']);

fclose(fp_debug);
end
