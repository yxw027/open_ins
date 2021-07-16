%show_result openrtk offline solution
%input:日期文件夹
clear;
clc;
filesname = 'E:/data/INS_ODO/Process/0826/test3/';
filesnamelist = dir(filesname);
filesize = size(filesnamelist,1);
n_res = 0; % number of solutions
for i = 1 : (filesize - 2)
    %当前文件夹名称
    folder  = [filesname filesnamelist(i+2).name '/'];
    filesnamelist1 = dir(folder);
    filesize1 = size(filesnamelist1,1);
    for j = 1 : (filesize1 - 2)
        if(strfind(filesnamelist1(j+2).name,'ref'))
            find_ref =1;
            refnamefolder = [folder filesnamelist1(j+2).name '/'];
            refnamefolderlist = dir(refnamefolder);
            cfg.refname = [refnamefolder refnamefolderlist(3).name];
            break;
        end
    end
    for j = 1 : (filesize1 - 2)
        if(strfind(filesnamelist1(j+2).name,'time.txt'))
            cfg.cifgname = [folder filesnamelist1(j+2).name];
            [ret,cfg.cfg_info] = readconfig(cfg.cifgname);
            if (1 == ret)
            find_config =1;
            cfg.cfg_info.legend = {'a', 'b', 'c'};
            end         
            break;
        end
    end
    

    if(1 == find_ref && 1 == find_config)
        for j = 1 : (filesize1 - 2)
            if(strfind(filesnamelist1(j+2).name,'dev'))
                devnamefolder = [folder filesnamelist1(j+2).name '/'];
                devnamefolderlist = dir(devnamefolder);
                filesize2 = size(devnamefolderlist,1);
                 need_to_be =  0;
                for k = 1 : (filesize2 - 2)
                    %需要分析的单次测试数据
                     if(strfind(devnamefolderlist(k+2).name,'-gnssposvel.txt'))
                         cfg.gnssfilename = [devnamefolder devnamefolderlist(k+2).name];
                         need_to_be =need_to_be + 1;
                         continue;
                     end
                     if(strfind(devnamefolderlist(k+2).name,'-imu.txt'))
                         cfg.imufilename = [devnamefolder devnamefolderlist(k+2).name];
                         need_to_be =need_to_be + 1;
                         continue;
                     end
                     if(strfind(devnamefolderlist(k+2).name,'_result.txt'))
                         cfg.solfilename = [devnamefolder devnamefolderlist(k+2).name];
                         need_to_be =need_to_be + 1;
                         continue;
                     end
                     if(strfind(devnamefolderlist(k+2).name,'_deb.txt'))
                         cfg.filfilename = [devnamefolder devnamefolderlist(k+2).name];
                         need_to_be =need_to_be + 1;
                         continue;
                     end
                end
                if (4 == need_to_be)
                    n_res = n_res + 1;
                    %填入需要输入的文件名，以及需要输出的文件夹
                    k = strfind(filesname,'/');
                    cfg.outputfilefolder = [filesname(1:k(end-1)) 'result' filesname(k(end-1):k(end)) filesnamelist(i+2).name  '/' filesnamelist1(j+2).name '/'];
                    show_result_all(n_res, cfg.imufilename,cfg.gnssfilename,cfg.solfilename,cfg.filfilename,cfg.outputfilefolder,cfg.refname,cfg.cfg_info); %,cfg.cfg_info
                end          
            end
        end
    else
        break;
    end
    
    
end
