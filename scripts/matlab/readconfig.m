function [ret,cfg_info] = readconfig(cifgname)
fid_conf = fopen(cifgname,'rt');
if fid_conf == -1
    disp('Cannont open config file');
    ret = -1;
    return;
end
cout = 0;
cfg_info.cdf_fh = [];
cfg_info.cep_fh = [];

while ~feof(fid_conf)
    line = fgetl(fid_conf);
    if(contains(line,'ref2'))
        if(contains(line,'gnss'))
            start_n = strfind(line,':');
            end_n = strfind(line,';');
            mid_n = strfind(line,',');
            cfg_info.leverarm.ref2gnss = [str2double(line(start_n(1)+1:mid_n(1)-1));...
                str2double(line(mid_n(1)+1:mid_n(2)-1));str2double(line(mid_n(2)+1:end_n(1)-1))];
            continue;
        end
        if(contains(line,'gnss'))
            start_n = strfind(line,':');
            end_n = strfind(line,';');
            mid_n = strfind(line,',');
            cfg_info.leverarm.ref2gnss = [str2double(line(start_n(1)+1:mid_n(1)-1));...
                str2double(line(mid_n(1)+1:mid_n(2)-1));str2double(line(mid_n(2)+1:end_n(1)-1))];
            continue;
        end
        if(contains(line,'imu'))
            start_n = strfind(line,':');
            end_n = strfind(line,';');
            mid_n = strfind(line,',');
            cfg_info.leverarm.ref2imu = [str2double(line(start_n(1)+1:mid_n(1)-1));...
                str2double(line(mid_n(1)+1:mid_n(2)-1));str2double(line(mid_n(2)+1:end_n(1)-1))];
            continue;
        end
        if(contains(line,'ublox'))
            start_n = strfind(line,':');
            end_n = strfind(line,';');
            mid_n = strfind(line,',');
            cfg_info.leverarm.ref2imu = [str2double(line(start_n(1)+1:mid_n(1)-1));...
                str2double(line(mid_n(1)+1:mid_n(2)-1));str2double(line(mid_n(2)+1:end_n(1)-1))];
            continue;
        end
        continue;
    end
    if(contains(line,':'))
        clear start_n end_n mid_n;
        cout = cout + 1;
        start_n = strfind(line,':');
        end_n = strfind(line,';');
        mid_n = strfind(line,',');
        cfg_info.dsc_nogs(cout) = string(line(2:start_n(1)-2));
%         cfg_info.cdf_fh(cout) = figure(100+cout);
%         cfg_info.cep_fh(cout) = figure(200+cout);
        cfg_info.t_nogps(cout,1) = str2double(line(start_n(1)+1:mid_n(1)-1));
        cfg_info.t_nogps(cout,2) = str2double(line(mid_n(1)+1:end_n(1)-1));
        continue;
    end
end
ret = 1;
end