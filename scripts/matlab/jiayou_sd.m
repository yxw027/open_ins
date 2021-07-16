% analyze OpenRTK330LI post-processing solution
% input: 
clear;
% clc;
close all
filesname = 'E:/data/INS_ODO/Process/0827/test1/';
filesnamelist = dir(filesname);
filesize = size(filesnamelist,1);
find_ref = 0; find_config = 0;

n_res = 0; % number of solutions

for i = 1 : (filesize - 2)
    filesnamelist(i+2).name
    if(strfind(filesnamelist(i+2).name,'ref'))
        find_ref =1;
        refnamefolder = [filesname filesnamelist(i+2).name '/'];
        refnamefolderlist = dir(refnamefolder);
        cfg.refname = [refnamefolder refnamefolderlist(3).name]; 

    elseif(strfind(filesnamelist(i+2).name,'time.txt'))
        cfg.cifgname = [filesnamelist(i+2).name];
        [ret,cfg.cfg_info] = readconfig([filesname cfg.cifgname]);
        if (1 == ret)
            find_config =1;
            cfg.cfg_info.legend = {'a'};
        end         
    end
    
    if(1 == find_ref && 1 == find_config)
       break; 
    end
    
end
    
if(1 == find_ref && 1 == find_config)
   for i = 1 : (filesize - 2)
      if(strfind(filesnamelist(i+2).name,'dev'))
        devnamefolder = [filesname filesnamelist(i+2).name '/'];
        devnamefolderlist = dir(devnamefolder);
        filesize2 = size(devnamefolderlist,1);
        need_to_be = 0;
        for k = 1 : (filesize2 - 2)
                    
             if(strfind(devnamefolderlist(k+2).name,'-gnssposvel.txt'))
                 cfg.gnssfilename = [devnamefolder devnamefolderlist(k+2).name];
                 need_to_be =need_to_be + 1;
                 continue;
             end
             if(strfind(devnamefolderlist(k+2).name,'-imu.txt'))
                 cfg.imufilename = [devnamefolder devnamefolderlist(k+2).name];
                 need_to_be =need_to_be + 1;
                 continue;
             else
                 cfg.imufilename = [];
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
         if (need_to_be >= 3)
            n_res = n_res + 1;
            k = strfind(filesname,'/');
            cfg.outputfilefolder = [filesname(1:k(end-1)) filesname(k(end-1):k(end)) 'result/' filesnamelist(i+2).name  '/'];
            show_result_all(n_res, cfg.imufilename,cfg.gnssfilename,cfg.solfilename,cfg.outputfilefolder,[], cfg.refname,cfg.cfg_info); %,cfg.cfg_info
         else
             disp('Missing input files');
         end          
    end
   end
else
    disp('No ref or config file found');
end

    
% end
