% Modified by Brandy Armstrong January 2012 to use only NCTOOLBOX 
% and Matlab builtin functions to read and write netcdf files
% jcw Feb 2019 - only use matalb BI
% modified by pzy, Sep 27, 2023

% 该函数是用于
% 该函数有以下输入参数：
% 1.T1：生成clm数据的起始日期
% 2.gn：一个结构体数组，里面包含ROMS自定义网格中的所有数据。详见get_roms_grid
% 3.clm：包含hycom索引的结构体
% 4.wdr：当前工作目录
% 5.clmname = climatology文件的网格名前缀
% 6.url：数据下载地址

function [fn]=updatclim_coawst_mw(T1, gn, clm, clmname, wdr, url)

% 首先，需要确定插值时间段的索引

disp('正在获取时间记录数...');

%{
以下四行代码不可用，需做出适当修改。最新的HYCOM数据集基准时间为2000年1月1日
00:00:00，时间分辨率为3小时，并且没有变量"MT"。

t0是2000年1月1日的MATLAB日期数字。time是从URL下载的时间序列数组，表示从2018年12
月4日12时至2023年9月30日21时的数据。由于这个时间序列是相对于2000年1月1日的，所以
我们需要用tg = time/24 + t0将下载的时间序列转换为MATLAB日期数字，即相对于公元零年
（MATLAB 的默认日期基准）的天数。这就是tg的含义。注意，因为数据集的时间分辨率为每
3小时，因此必须要先将time除以24再和时间基准t0相加，否则会得到错误的结果。

然后，使用datestr函数将tg转换为字符串形式的日期和时间，然后再用 str2num 函数将其
转换为数字。这些数字被用作 julian 函数的输入参数，该函数返回儒略日期（即相对于公元
前4713年1月1日的天数）。这就是tg2的含义。使用天文惯例将公历日期转换为十进制儒略日，
但时间零从午夜开始，而不是从正午开始。在这个惯例中，朱利安日2440000从1968年5月23日
0时开始。十进制儒略历日，与Matlab的双精度，产生的十进制日的精度约为0.1毫秒。
%}
% t0 = datenum(1900,12,31); % tr0=datenum(1858,11,17);
% time = ncread(url,'MT');
% tg=time+t0;
% tg2=julian(str2num(datestr(tg,'yyyy')),str2num(datestr(tg,'mm')),str2num(datestr(tg,'dd')),str2num(datestr(tg,'HH')))-2400001;
t0 = datenum(2000,01,01,0,0,0);
time = ncread(url,'time',1,14085);    % 下载从2018年12月4日12时至2023年9月30日21时的数据
tg = time/24 + t0;
tg2 = julian(str2num(datestr(tg,'yyyy')),str2num(datestr(tg,'mm')),...
    str2num(datestr(tg,'dd')),str2num(datestr(tg,'HH')))-2400001;


% 获取用户定义的初始时间

%{
接下来四行代码的逻辑是：寻找hycom数据集时间tg和用户定义起始时间T1的交集。intersect
函数的作用是返回两个数组的公共元素，以及它们在各自数组中的索引。在这里，junk是公共
元素，tid1是公共元素在tg中的索引，ib是公共元素在floor(T1)中的索引。而floor(T1)
是将T1中的每个元素向下取整到最接近的整数。

if isempty(tid1)是为了检查 tid1 是否为空。如果intersect函数没有找到任何公共元素
，那么就将tid1设置为tg的长度,也即tg最后一个元素的索引。这是为了在后续的代码中使用
tid1 作为索引，确保其总是有效的而不会报错。该段近一次被jc warner于2014年8月23日
修改过。
%}
% [junk,tid1,ib]=intersect(tg,floor(T1));
[~,tid1,~]=intersect(tg,floor(T1));
if isempty(tid1)
  tid1=length(tg);
end

%{
这里将输入的clmname(即climatology文件的网格名前缀)赋给fn，然后调用
create_roms_netcdf_clm_mwUL生成对应的nc文件。
%}
fn = clmname;
disp(['正在创建nc文件',fn,'...']);
create_roms_netcdf_clm_mwUL(fn,gn,1);% converted to BI functions

%使用builtin (BI) functions填充网格维度

%{
这里使用ncid = netcdf.open(source,mode)以mode指定的访问类型打开source所指向的nc
文件。mode的值可取'WRITE'、'SHARE' 或 'NOWRITE'。

然后，使用varid = netcdf.inqVarID(ncid,varname)返回与varname对应的varid。varname
指定变量名称。ncid指定netcdf.create或netcdf.open返回的netCDF文件标识符。

详见netcdf.open和netcdf.inqVarID的说明文档。
%}
RN = netcdf.open(fn,'WRITE');
lonid = netcdf.inqVarID(RN,'lon_rho');
netcdf.putVar(RN,lonid,gn.lon_rho);
latid = netcdf.inqVarID(RN,'lat_rho');
netcdf.putVar(RN,latid,gn.lat_rho);
netcdf.close(RN)

%%

%{
repmat函数用于将某个矩阵按照行方向和列方向重复若干次，从而生成一个新的矩阵。例如，
在这里clm.lon是一个102*1的数组，clm.lat是一个171*1的数组；那么，
X = repmat(clm.lon,1,length(clm.lat))表示将矩阵clm.lon在行方向上的重复次数*1
(即无需重复)，在列的方向上重复次数*171(即重复170次)，这样得到的数组X的大小是102*171。
%}
tz_levs = length(clm.z);
X = repmat(clm.lon,1,length(clm.lat));
Y = repmat(clm.lat,length(clm.lon),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['正在为',datestr(tg(tid1)),'处的clm文件插值water_u...']);

%{
下面开始调用插值函数进行插值。为clm结构体新设一个变量u，其大小为40(hycom垂直层级数)
*836(ROMS自定义网格的L)*753(ROMS自定义网格的M)，
%}
ttu=1;
clm.u=zeros([length(clm.z) size(gn.lon_rho)]);
while ttu==1
    try
        tmpt=ncread(url,'water_u',[clm.ig0 clm.jg0 1 tid1],[clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        for k=1:tz_levs
            disp(['doing griddata u for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            F = scatteredInterpolant(X(:),Y(:),tmp(:));
            cff = F(gn.lon_rho,gn.lat_rho);
            clm.u(k,:,:)=maplev(cff);
        end
        ttu=0;
    catch
        disp(['catch u Unable to download HYCOM u data at' datestr(now)]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM u data at');
        fprintf(fid,datestr(now));
        fprintf(fid,'\n');
    end
end

% 实现从标准z坐标到s坐标的垂向插值(t,s,u,v)

u = roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.u,gn,'u',0);
clm = rmfield(clm,'u');
save u.mat u
clear u;

disp(['正在为',datestr(tg(tid1)),'处的clm文件插值water_v...']);
ttv=1;
clm.v=zeros([length(clm.z) size(gn.lon_rho)]);
while ttv==1
    try
        tmpt=ncread(url,'water_v',[clm.ig0 clm.jg0 1 tid1],[clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        for k=1:tz_levs
            disp(['doing griddata v for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            F = scatteredInterpolant(X(:),Y(:),tmp(:));
            cff = F(gn.lon_rho,gn.lat_rho);
            clm.v(k,:,:)=maplev(cff);
        end
        ttv=0;
    catch
        disp(['catch v Unable to download HYCOM v data at' datestr(now)]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM v data at');
        fprintf(fid,datestr(now));
        fprintf(fid,'\n');
    end
end
%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
v = roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.v,gn,'v',0);
clm=rmfield(clm,'v');
save v.mat v
clear v;

%== Rotate the velocity
theta=exp(-sqrt(-1)*mean(mean(gn.angle)));
eval(['load ',wdr,'\','u.mat']);
eval(['load ',wdr,'\','v.mat']);
disp('doing rotation to grid for u and v');
uv=(u2rho_3d_mw(u)+sqrt(-1)*v2rho_3d_mw(v)).*theta;
u=rho2u_3d_mw(real(uv)); v=rho2v_3d_mw(imag(uv));
clear uv

%% == output
RN=netcdf.open(fn,'WRITE');

tempid=netcdf.inqVarID(RN,'u');
netcdf.putVar(RN,tempid,shiftdim(u,1));

tempid=netcdf.inqVarID(RN,'v');
netcdf.putVar(RN,tempid,shiftdim(v,1));

clear u; clear v;
tempid=netcdf.inqVarID(RN,'ocean_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'zeta_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'v2d_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'v3d_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'salt_time');
netcdf.putVar(RN,tempid,tg2(tid1));
tempid=netcdf.inqVarID(RN,'temp_time');
netcdf.putVar(RN,tempid,tg2(tid1));
netcdf.close(RN);
%%
%== Depth averaging u, v to get Ubar
eval(['load ',wdr,'\','u.mat']);
eval(['load ',wdr,'\','v.mat']);
cc=roms_zint_mw(u,gn);  ubar=rho2u_2d_mw(u2rho_2d_mw(cc)./gn.h);
cc=roms_zint_mw(v,gn);  vbar=rho2v_2d_mw(v2rho_2d_mw(cc)./gn.h);
%== Rotate the velocity
uv=(u2rho_2d_mw(ubar)+sqrt(-1)*v2rho_2d_mw(vbar)).*theta;
ubar=rho2u_2d_mw(real(uv)); vbar=rho2v_2d_mw(imag(uv));
clear u
clear v

RN=netcdf.open(fn,'WRITE');
tempid=netcdf.inqVarID(RN,'ubar');
netcdf.putVar(RN,tempid,ubar);
tempid=netcdf.inqVarID(RN,'vbar');
netcdf.putVar(RN,tempid,vbar);
netcdf.close(RN);

clear ubar
clear vbar
clear uv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the zeta data
disp(['正在为',datestr(tg(tid1)),'处的clm文件插值zeta(自由表面)...']);
ttz=1;
while ttz==1
    try
        tmpt=ncread(url,'surf_el',[clm.ig0 clm.jg0 tid1],[clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 1 ] );
        tmp=double(squeeze(tmpt(:,:)));
        disp('doing griddata zeta for HYCOM ');
        F = scatteredInterpolant(X(:),Y(:),tmp(:));
        cff = F(gn.lon_rho,gn.lat_rho);
        zeta=maplev(cff);
        ttz=0;
    catch
        disp(['catch z Unable to download HYCOM ssh data at' datestr(now)]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM ssh data at');
        fprintf(fid,datestr(now));
        fprintf(fid,'\n');
    end
end
clear tmp
%
%== output zeta
%
RN=netcdf.open(fn,'WRITE');
tempid=netcdf.inqVarID(RN,'zeta');
netcdf.putVar(RN,tempid,zeta);
netcdf.close(RN);
clear zeta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['正在为',datestr(tg(tid1)),'处的clm文件插值温度...']);
ttt=1;
clm.temp=zeros([length(clm.z) size(gn.lon_rho)]);
while ttt==1
    try
        tmpt=ncread(url,'water_temp',[clm.ig0 clm.jg0 1 tid1],[clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        for k=1:tz_levs
            disp(['doing griddata temp for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            F = scatteredInterpolant(X(:),Y(:),tmp(:));
            cff = F(gn.lon_rho,gn.lat_rho);
%           cff(cff<0)=nan;
            clm.temp(k,:,:)=maplev(cff);
        end
        ttt=0;
    catch
        disp(['catch temp Unable to download HYCOM temp data at' datestr(now)]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM temp data at');
        fprintf(fid,datestr(now));
        fprintf(fid,'\n');
    end
end
%
%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
%
temp=roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.temp,gn,'rho',0);
clm=rmfield(clm,'temp');
%
%== output temp
%
RN=netcdf.open(fn,'WRITE');
tempid=netcdf.inqVarID(RN,'temp');
netcdf.putVar(RN,tempid,shiftdim(temp,1));
netcdf.close(RN);
clear temp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['正在为',datestr(tg(tid1)),'处的文件插值盐度...']);
tts=1;
clm.salt=zeros([length(clm.z) size(gn.lon_rho)]);
while tts==1
    try
        tmpt=ncread(url,'salinity',[clm.ig0 clm.jg0 1 tid1],[clm.ig1-clm.ig0+1 clm.jg1-clm.jg0+1 tz_levs 1 ] );
        for k=1:tz_levs
            disp(['doing griddata salt for HYCOM level ' num2str(k)]);
            tmp=double(squeeze(tmpt(:,:,k)));
            F = scatteredInterpolant(X(:),Y(:),tmp(:));
            cff = F(gn.lon_rho,gn.lat_rho);
            cff(cff<0)=nan;
            clm.salt(k,:,:)=maplev(cff);
        end
        tts=0;
    catch
        disp(['catch temp Unable to download HYCOM temp data at' datestr(now)]);
        fid=fopen('coawstlog.txt','a');
        fprintf(fid,'Unable to download HYCOM temp data at');
        fprintf(fid,datestr(now));
        fprintf(fid,'\n');
    end
end
%
%== Vertical interpolation (t,s,u,v) from standard z-level to s-level
%
salt = roms_from_stdlev_mw(gn.lon_rho,gn.lat_rho,clm.z,clm.salt,gn,'rho',0);
clm = rmfield(clm,'salt');
%
%== output salt
%
RN=netcdf.open(fn,'WRITE');
tempid=netcdf.inqVarID(RN,'salt');
netcdf.putVar(RN,tempid,shiftdim(salt,1));
netcdf.close(RN);
clear salt;

disp(['于',datestr(now),'完成clm文件创建。']);
%%
