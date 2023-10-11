% roms_master_climatology_coawst_mw
%
% This routine :
%  - creates climatology, boundary, and initial condition files for ROMS: 
%    coawst_clm.nc ; coawst_bdy.nc ; coawst_ini.nc 
%    on a user-defined grid for a user-defined date.
%
% This is currently set up to use opendap calls to acquire data
% from HYCOM + NCODA Global 1/12 Degree Analysis and interp to roms grid.
%  
% based on efforts by:
% written by Mingkui Li, May 2008
% Modified by Brandy Armstrong March 2009
% jcwarner April 20, 2009
% Ilgar Safak modified on June 27, 2012 such that now:
% - HYCOM url is a user-definition
% - "hc" is called from the structure "gn".(still needs to be tested with wet/dry).
% - updatinit_coawst_mw.m modified to get desired time (T1) as a variable;
%    ocean_time=T1-datenum(1858,11,17,0,0,0)
% Updates from Christie Hegermiller, Feb 2019

% modified by pzy, Sep 27, 2023

% 该程序用于获取气候数据、初始场和边界的nc格式文件。使用的数据为HYCOM和NCODA的全球1/12°的
% 分析数据，HYCOM url由自己定义。hc参数由gn调用；可以更改updatinit_coawst_mw.m来获取合适
% 的T1值作为初始时间。

%%%%%%%%%%%%%%%%%%%%%   从用户输入开始  %%%%%%%%%%%%%%%%%%%%%%%%%%

% (1)首先，输入开始时间T1和需要获得clm数据的天数

T1 = datenum(2022,07,01,12,0,0);     % 定义模式开始运行的起始时间，这里是从2022年7月1日12时开始

% 这里代表用户想运行的天数，以及生成clm文件的频率

numdays = 5;                         % 定义模式运行的天数，这里为30天
dayFrequency = 1;                    % 定义每隔几天生成一次clm数据，这里设置为1天

% (2) 然后，键入HYCOM目录中对应开始时间T1所指向的下载URL地址。
% 可见http://tds.hycom.org/thredds/catalog.html
 
url = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0';  % 2018-12-04至今

% (3) 键入工作目录

wdr = 'E:\研究生课件\COAWST\ROMS学习\upwelling_project';

% (4) 键入ROMS网格名称和路径

modelgrid = 'upwelling_roms_grid.nc';

% (5) 输一些入网格垂直坐标参数。这些参数需要与rom设置和自己的网格保持一致。

theta_s     = 8.0;
theta_b     = 4.0;
Tcline      = 20.0;
N           = 30;
Vtransform  = 2;
Vstretching = 4;

%%%%%%%%%%%%%%%%%%%%%   用户输入结束   %%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1.生成ROMS启动时间的clm、bdy和int文件

eval(['cd ',wdr])     % cd到工作目录中

t1_start = tic;                   % 启动定时器，记录下当前的时间。

% 调用get_ijrg，以获取自定义的roms网格的HYCOM索引

disp('正在获取roms网格、hycom网格和重叠部分的索引...') 
[gn, clm] = get_ijrg(url, modelgrid, theta_s, theta_b, Tcline, N,...
    Vtransform, Vstretching);

% 调用updatclim_coawst_mw生成clm文件

disp('正在创建clm文件...')                           % 正在生成clm文件
fn = updatclim_coawst_mw(T1, gn, clm, 'coawst_clm.nc', wdr, url);

% 调用updatbdry_coawst_mw生成bdy文件

disp('正在创建bdy文件...')                           % 正在生成bdy文件
updatbdry_coawst_mw(fn, gn, 'coawst_bdy.nc', wdr)

% 调用updatinit_coawst_mw生成ini文件

disp('正在创建初始场文件...')                         % 正在生成int文件         
updatinit_coawst_mw(fn, gn, 'coawst_ini.nc', wdr, T1)

t1_end = toc(t1_start);     % 调用toc来停止计时器计算时间差，测量代码的执行时间,评估代码的效率。


%% 2.生成长期的clm文件

t2_start = tic;

%{
接下来的if语句逻辑是：若模拟时间大于1天，则继续生成后续天数的clm和bdy文件，
否则什么也不做。
%}

if numdays>1          
    disp('正在创建多日的clm和边界文件...')

    %{
    接下来的if语句的逻辑是：检测当前的操作系统。根据不同的操作系统，使用不
    同的复制指令执行相同的操作：将前面生成的clm和bdy文件复制并重命名为
    coawst_clm_(自己选择的初始时间).nc。ispc语句可以判断当前自己的操作系统
    ，如果是windows则返回True，否则返回false。copy和cp分别是windows和linux
    的复制指令，前面加上"!"是为了让matlab知道要将该文本传递给操作系统进行处理。
    %}

    if (ispc)
      eval(['!copy coawst_clm.nc coawst_clm_', datestr(T1,'yyyymmdd'),'.nc'])
      eval(['!copy coawst_bdy.nc coawst_bdy_', datestr(T1,'yyyymmdd'),'.nc'])
    else
      eval(['!cp coawst_clm.nc coawst_clm_', datestr(T1,'yyyymmdd'),'.nc'])
      eval(['!cp coawst_bdy.nc coawst_bdy_', datestr(T1,'yyyymmdd'),'.nc'])
    end

    %{
    接下来的for循环的逻辑是：根据上面设置的dayFrequency参数，每隔若干天生成一个
    bdy和clm文件，并对这些文件用相应的日期命名。
    %}
    for it = dayFrequency:dayFrequency:numdays-1
        fname = ['coawst_clm_', datestr(T1+it,'yyyymmdd'),'.nc'];
        fn = updatclim_coawst_mw(T1+it,gn,clm,fname,wdr,url);
        fname = ['coawst_bdy_', datestr(T1+it,'yyyymmdd'),'.nc'];
        updatbdry_coawst_mw(fn,gn,fname,wdr)
    end
    
    %{
    接下来两行的逻辑是：用dirsort.m文件定义的函数来将生成的nc文件按日期排列整齐
    ，并得到存放它们信息的结构体数组。关于dirsort的逻辑详见dirsort.m。
    %}
    Dclm = dirsort('coawst_clm_*.nc');   
    Dbdy = dirsort('coawst_bdy_*.nc');

    fout = 'merged_coawst_clm.nc';      % 定义merged clm文件名
    foutb = 'merged_coawst_bdy.nc';     % 定义merged bdy文件名

    %create netcdf files to merge climatology into
    create_roms_netcdf_clm_mwUL(fout,gn,length(Dclm));% converted to BI functions
    create_roms_netcdf_bndry_mwUL(foutb,gn,length(Dbdy));% converted to BI functions
    %% fill merged climatology files with data from each clm file
    % each file must contain only ONE time step
    %get variable names
    vinfo = ncinfo(fout);
    for nf = 1:length(Dclm)
        fin = Dclm(nf).name;
        for nv = 1:length({vinfo.Variables.Name})
            if length({vinfo.Variables(nv).Dimensions.Name}) == 4
                eval(['ncwrite(fout,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 1 1 nf]);']);
            elseif length({vinfo.Variables(nv).Dimensions.Name})==3
                eval(['ncwrite(fout,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 1 nf]);']);
            elseif length({vinfo.Variables(nv).Dimensions.Name})==2
                try
                    eval(['ncwrite(fout,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 nf]);']);
                catch
                    display([vinfo.Variables(nv).Name ' is a dimension and has already been written to the file.'])
                end
            elseif length({vinfo.Variables(nv).Dimensions.Name})==1
                try
                    eval(['ncwrite(fout,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[nf]);']);
                catch
                    display([vinfo.Variables(nv).Name ' is a dimension and has already been written to the file.'])
                end
            end
        end
    end
    
    vinfo=ncinfo(foutb);
    for nf=1:length(Dbdy)
        for nv=1:length({vinfo.Variables.Name})
            fin=Dbdy(nf).name;
            if length({vinfo.Variables(nv).Dimensions.Name})==4
                eval(['ncwrite(foutb,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 1 1 nf]);']);
            elseif length({vinfo.Variables(nv).Dimensions.Name})==3
                eval(['ncwrite(foutb,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 1 nf]);']);
            elseif length({vinfo.Variables(nv).Dimensions.Name})==2
                try
                    eval(['ncwrite(foutb,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[1 nf]);']);
                catch
                    display([vinfo.Variables(nv).Name ' is a dimension and has already been written to the file.'])
                end
                
            elseif length({vinfo.Variables(nv).Dimensions.Name})==1
                try
                    eval(['ncwrite(foutb,''',vinfo.Variables(nv).Name,''',ncread(fin,''',vinfo.Variables(nv).Name,'''),[nf]);']);
                catch
                    display([vinfo.Variables(nv).Name ' is a dimension and has already been written to the file.'])
                end
            end
        end
    end
end

t2_end = toc(t2_start);

