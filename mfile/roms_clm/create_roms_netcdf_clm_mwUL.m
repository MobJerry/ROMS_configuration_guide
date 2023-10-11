%{
该函数用于生成前缀名称为fn的空白nc文件，以便之后把变量写入该文件。

该函数的输入参数有：
1.fn：生成的空白nc文件的前缀名
2.gn：包含该nc文件对应的ROMS网格全部变量信息的结构体
3.t_clim：nc文件变量的时间维度。

modified by pzy, Oct 10, 2023.
%}
function create_roms_netcdf_clm_mwUL(fn,gn,t_clim)

[xi,eta]=size(gn.lon_rho);        % 根据变量lon_rho的大小确定rho_point处的xi坐标(x方向)
                                  % 和eta坐标(y方向)的大小

%{
ncid = netcdf.create(filename,cmode)根据文件创建模式"cmode"新建netCDF文件。返回
值ncid是一个文件ID。

cmode 参数确定文件访问的类型，详见netcdf.create的matlab文档；而'CLOBBER'的含义是
覆盖现有的同名文件。

if isempty(nc), return, end 这行代码检查nc是否为空。如果nc为空(说明NetCDF文件创
建失败)，那么函数就会立即返回，不再执行后面的代码。
%}
nc = netcdf.create(fn,'clobber');
if isempty(nc), return, end

disp(' ## 正在定义全局变量...')

%{
netcdf.putAtt(ncid,varid,attrname,attrvalue)将名为"attrname"的属性及其值attrvalue
写入"varid"指定的netCDF变量。要指定全局属性，可令varid=netcdf.getConstant('NC_GLOBAL')。

ncid 是 netcdf.create 或 netcdf.open 返回的 NetCDF 文件标识符。

netcdf.putAtt(ncid,varid,attrname,attrvalue,xtype)将attrvalue作为在"xtype"中
指定的数据类型写入。"xtype"可指定为NC_DOUBLE、NC_FLOAT、NC_INT64等，详见netcdf.putAtt
的帮助文档。
%}
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', ...
    ['Created by updatclim on ' datestr(now)]);
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type', ...
    'climate forcing file from http://hycom.coaps.fsu.edu:8080/thredds/dodsC/glb_analysis');

% 接下来创建nc文件变量的维度

disp(' ## 正在定义数据维度...')

%{
接下来，使用dimid = netcdf.defDim(ncid,dimname,dimlen)在"ncid"指定的netCDF文件
中创建新维度。其中"dimname"是指定维度名称的字符向量或字符串标量，dimlen是指定其长
度的数值。要定义无限制维度，则令dimlen=netcdf.getConstant('NC_UNLIMITED')。

netcdf.defDim 返回与新维度相对应的数值ID dimid，代表该维度是第几个生成的。

详见netcdf.defDim的帮助文档。
%}
LP = xi;
MP = eta;
L = LP-1;
M = MP-1;
N = 0;

psidimID = netcdf.defDim(nc,'xpsi',L);
xrhodimID = netcdf.defDim(nc,'xrho',LP);
xudimID = netcdf.defDim(nc,'xu',L);
xvdimID = netcdf.defDim(nc,'xv',LP);

epsidimID = netcdf.defDim(nc,'epsi',M);
erhodimID = netcdf.defDim(nc,'erho',MP);
eudimID = netcdf.defDim(nc,'eu',MP);
evdimID = netcdf.defDim(nc,'ev',M);
s_rhodimID = netcdf.defDim(nc,'s_rho',gn.N);

octdimID = netcdf.defDim(nc,'ocean_time',t_clim);
zttdimID = netcdf.defDim(nc,'zeta_time',t_clim);
v2tdimID = netcdf.defDim(nc,'v2d_time',t_clim);
v3tdimID = netcdf.defDim(nc,'v3d_time',t_clim);
sltdimID = netcdf.defDim(nc,'salt_time',t_clim);
tptdimID = netcdf.defDim(nc,'temp_time',t_clim);
onedimID = netcdf.defDim(nc,'one',1);

% 接下来，写入变量和属性Variables and attributes:

disp(' ## 正在定义变量和属性...')

%{
varid = netcdf.defVar(ncid,varname,xtype,dimids)在ncid所标识的nc文件中创建一个
新变量。"varname"指定变量名称，"xtype"指定变量的NetCDF数据类型。"dimids"指定前面
已经创建的维度的ID列表。

netcdf.defVar返回varid，新变量的数值标识符。

详见的netcdf.defVar帮助文档。
%netcdf.putAtt(ncid,varid,attrname,attrvalue)
%}
ocID = netcdf.defVar(nc,'ocean_time','double',octdimID);
netcdf.putAtt(nc,ocID,'long_name','wind field time');
netcdf.putAtt(nc,ocID,'units','days');
netcdf.putAtt(nc,ocID,'field','wave_time, scalar, series');

ztID = netcdf.defVar(nc,'zeta_time','double',zttdimID);
netcdf.putAtt(nc,ztID,'long_name','zeta_time');
netcdf.putAtt(nc,ztID,'units','days');
netcdf.putAtt(nc,ztID,'field','zeta_time, scalar, series');

v2ID = netcdf.defVar(nc,'v2d_time','double',v2tdimID);
netcdf.putAtt(nc,v2ID,'long_name','v2d_time');
netcdf.putAtt(nc,v2ID,'units','days');
netcdf.putAtt(nc,v2ID,'field','v2d_time, scalar, series');

v3ID = netcdf.defVar(nc,'v3d_time','double',v3tdimID);
netcdf.putAtt(nc,v3ID,'long_name','v3d_time');
netcdf.putAtt(nc,v3ID,'units','days');
netcdf.putAtt(nc,v3ID,'field','v3d_time, scalar, series');

slID = netcdf.defVar(nc,'salt_time','double',sltdimID);
netcdf.putAtt(nc,slID,'long_name','salt_time');
netcdf.putAtt(nc,slID,'units','days');
netcdf.putAtt(nc,slID,'field','salt_time, scalar, series');

tpID = netcdf.defVar(nc,'temp_time','double',tptdimID);
netcdf.putAtt(nc,tpID,'long_name','temp_time');
netcdf.putAtt(nc,tpID,'units','days');
netcdf.putAtt(nc,tpID,'field','temp_time, scalar, series');

lonID = netcdf.defVar(nc,'lon_rho','float',[xrhodimID erhodimID]);
netcdf.putAtt(nc,lonID,'long_name','lon_rho');
netcdf.putAtt(nc,lonID,'units','degrees');
netcdf.putAtt(nc,lonID,'FillValue_',100000.);
netcdf.putAtt(nc,lonID,'missing_value',100000.);
netcdf.putAtt(nc,lonID,'field','xp, scalar, series');

latID = netcdf.defVar(nc,'lat_rho','float',[xrhodimID erhodimID]);
netcdf.putAtt(nc,latID,'long_name','lon_rho');
netcdf.putAtt(nc,latID,'units','degrees');
netcdf.putAtt(nc,latID,'FillValue_',100000.);
netcdf.putAtt(nc,latID,'missing_value',100000.);
netcdf.putAtt(nc,latID,'field','yp, scalar, series');

zetID = netcdf.defVar(nc,'zeta','double',[xrhodimID erhodimID zttdimID]);
netcdf.putAtt(nc,zetID,'long_name','zeta');
netcdf.putAtt(nc,zetID,'units','meter');
netcdf.putAtt(nc,zetID,'field','zeta, scalar, series');

salID = netcdf.defVar(nc,'salt','float',[xrhodimID erhodimID s_rhodimID sltdimID]);
netcdf.putAtt(nc,salID,'long_name','salt');
netcdf.putAtt(nc,salID,'units','psu');
netcdf.putAtt(nc,salID,'field','salt, scalar, series');

tmpID = netcdf.defVar(nc,'temp','float',[xrhodimID erhodimID s_rhodimID tptdimID]);
netcdf.putAtt(nc,tmpID,'long_name','temp');
netcdf.putAtt(nc,tmpID,'units','C');
netcdf.putAtt(nc,tmpID,'field','temp, scalar, series');

uID = netcdf.defVar(nc,'u','float',[xudimID eudimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc,uID,'long_name','velx');
netcdf.putAtt(nc,uID,'units','meter second-1');
netcdf.putAtt(nc,uID,'field','velx, scalar, series');

vID = netcdf.defVar(nc,'v','float',[xvdimID evdimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc,vID,'long_name','vely');
netcdf.putAtt(nc,vID,'units','meter second-1');
netcdf.putAtt(nc,vID,'field','vely, scalar, series');

ubID = netcdf.defVar(nc,'ubar','float',[xudimID eudimID v2tdimID]);
netcdf.putAtt(nc,ubID,'long_name','mean velx');
netcdf.putAtt(nc,ubID,'units','meter second-1');
netcdf.putAtt(nc,ubID,'field','mean velx, scalar, series');

vbID = netcdf.defVar(nc,'vbar','float',[xvdimID evdimID v2tdimID]);
netcdf.putAtt(nc,vbID,'long_name','mean vely');
netcdf.putAtt(nc,vbID,'units','meter second-1');
netcdf.putAtt(nc,vbID,'field','mean vely, scalar, series');

netcdf.close(nc)

