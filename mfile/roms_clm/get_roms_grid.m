%
% 本m文件get_roms_grid的目的是构建或更新一个包含自定义ROMS网格所有变量信息的一个
% 结构体数组。 
%
% 该函数的用法有：
% 
% Gout = get_roms_grid(Ginp)
% Gout = get_roms_grid(Ginp, Sinp)
% Gout = get_roms_grid(Ginp, Sinp, Tindex)
% Gout = get_roms_grid(Ginp, Tindex)
%
% 该函数的输入参数详解如下:
%
% Ginp：可以输入以下两种形式：1.包含所有自定义网格变量的ROMS网格/历史nc文件/URL
% 的名称，以字符串的形式输入。例如，本地有一个ROMS网格nc文件，则可以将该网格的路径
% 名作为输入；2.一个包含垂直坐标的网格结构数组。
%
% Sinp：可接受以下两种输入：1.ROMS输出的nc文件/URL名称，从中可以确定出垂直坐标，          
% 以字符串的形式输入。2.一个结构数组，要包含以下垂直坐标的拉伸参数：
% Sinp.N               垂直层级数
% Sinp.Vtransform      垂直转换算法
% Sinp.Vstretching     垂直拉伸算法
% Sinp.theta_s         表层控制参数
% Sinp.theta_b         底层控制参数
% Sinp.Tcline          人为定义的Tcline
% Sinp.hc              ROMS所用的拉伸控制深度
%
% Tindex：时间记录指数，用于获取垂直坐标的自由表面信息。该参数是可以缺省的，输入为
% 一个标量。若Tindex缺省、为0或为空，会自动假设为zeta=0处没有扰动；否则，处理输入
% nc文件中的自由表面记录。自由表面可以从Ginp或Sinp(如果提供历史nc文件)中读取。
% 
% 该函数的输出参数如下：
% Gout：包含所有ROMS网格变量的结构数组。
%
% 例子:
%
% Gout = get_roms_grid('ocean_grd.nc');
%
% 会返回除了那些与垂直网格有关的变量之外的所有网格变量，以结构数组的形式保存。垂直
% 网格参数和自由表面的信息不在其中。
%
% 
% Gout = get_roms_grid('ocean_grd.nc', 'ocean_his.nc');
%
% 会返回包含所有网格变量的结构数组。在这种情况下垂向深度是没有扰动的，因为没有
% 指定Tindex。
%
%
% Gout = get_roms_grid('ocean_grd.nc', 'ocean_his.nc', MyTimeRec);
%
% 返回一个包含所有网格变量的结构数组。垂直深度与时间有关，因为指定了Tindex。
%
%
% Gout = get_roms_grid('ocean_his.nc', MyTimeRec);
%
% 返回一个包含所有网格变量的结构数组。这里作者假设ROMS历史文件"ocean_his.nc"
% 包含了所有的网格变量，所以没有激活NO_WRITE_GRID选项。
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
% modified by pzy,Sep 27, 2023

% 一.检查输入参数。

function Gout = get_roms_grid(Ginp, Sinp, Tindex)
  
process.horizontal = false;   % 定义结构数组process，给一些变量默认赋值为false
process.parameters = false;
process.zeta       = false;
process.vertical   = false;

TimeRecord = [];              % 定义空数组TimeRecord

is3d = false;                 % 定义一个判断网格是否三维的变量，默认赋值为false

%{ 
以下if语句的逻辑是：处理输入变量Ginp，并根据其类型（结构体或非结构体）执行不同的操作。

isstruct()是一个函数，用于检查输入是否为结构体。如果输入是结构体，它会返回1（真），否
则返回0（假）。而"~"是一个逻辑运算符，用于对逻辑值进行取反操作。也就是说，如果
isstruct(Ginp)返回1，那么~isstruct(Ginp)就会返回0，反之亦然。

所以，~isstruct(Ginp)的作用就是检查Ginp是否不是一个结构体。如果Ginp不是结构体，那么
~isstruct(Ginp)将返回1（真），否则返回0（假）。

若Ginp不是一个结构体数组，则：
1.将Gout的grid_name赋值为Ginp所指的名称；
2.使用nc_inq(Ginp)获取关于Ginp的信息，并将结果存储在Ginfo中；
3.从Ginfo.Variables.Name中提取所有变量名，并将其存储在vnames中；
4.从Ginfo.Attributes.Name中提取所有全局属性名，并将其存储在anames中；
5.计算变量的数量并将其存储在NvarsG中；
6.检查vnames中是否存在名为s_rho的变量。若存在，则将其视为三维ROMS网格应用，并将结
果存储在is3d中。

strcmp是MATLAB中的一个函数，用于比较字符串。这个函数会返回一个逻辑数组，数组中的每
个元素对应于vnames中的一个元素。如果vnames中的元素与s_rho相等，那么对应的数组元素
就为1（真），否则为0（假）。

而any用于检查输入数组中是否存在任何非零元素（在逻辑上等同于"真"）。如果strcmp(vnames,'s_rho')
返回的数组中存在任何1，那么any(strcmp(vnames,'s_rho'))就会返回1（真），否则返回
0（假）。所以，这行代码的作用就是检查变量名列表vnames中是否存在名为s_rho的变量。
如果存在，那么它将返回1（真），否则返回0（假）。

若Ginp是一个结构数组，则通过~isempty(Ginp.N)检查其是否包含字段垂直层级N。如果包含，
那么网格被视为3D应用。
%}
if (~isstruct(Ginp))
  Gout.grid_name = Ginp;
  Ginfo  = nc_inq(Ginp);
  vnames = {Ginfo.Variables.Name};
  anames = {Ginfo.Attributes.Name};          
  NvarsG = length(Ginfo.Variables);          
  is3d   = any(strcmp(vnames,'s_rho'));      
else
  is3d   = ~isempty(Ginp.N);
end

%{
接下来if语句的逻辑是：处理输入参数Sinp和Tindex，并根据其类型和数量执行不同的操作。

如果函数接收到超过一个输入参数（由nargin > 1判断），代码将执行以下操作：

1.如果Sinp是字符型（由ischar(Sinp)判断），那么将process.zeta设为真，并将Sinp赋
值给Gout.his_name。这时，Sinp被视为一个NetCDF文件。

2，如果Sinp是结构体（由isstruct(Sinp)判断），那么将process.parameters设为真。这里，
Sinp被视为一个包含垂直参数的结构体。

3.如果Sinp是一个数值，且函数只接收到两个输入参数（由isnumeric(Sinp) && nargin == 2判断）
，那么将process.zeta设为真，并将Sinp赋值给TimeRecord。这里，Sinp被视为Tindex。
%}
if (nargin > 1)
  if (ischar(Sinp))
    process.zeta = true;
    Gout.his_name = Sinp;
  elseif (isstruct(Sinp))
    process.parameters = true;               
  elseif (isnumeric(Sinp) && nargin == 2)
    process.zeta = true;
    TimeRecord = Sinp;                       
  end
end

%{
如果函数接收到三个输入参数（由 nargin == 3 判断），代码将执行以下操作：将 Tindex
赋值给 TimeRecord。
%}
if (nargin == 3)
  TimeRecord = Tindex;
end

% 二.初始化

spherical = false;           % 定义一个判断网格是否在球面上的变量，默认赋值为false

%{
以下这段很长的嵌套if语句的主要目的是判断Ginp的类型，并根据不同类型执行不同的操作。

如果Ginp不是结构体（由~isstruct(Ginp)判断），代码将执行以下操作：

1.检查全局属性名列表anames中是否存在名为coarse_factor或refine_factor的属性。如果
存在，那么检查是否存在名为parent_grid的属性，并将其值赋给Gout.parent_grid；如果不
存在，什么也不做。

2.将一些字段（如TimeRecord, Lm, Mm, N, coarse_factor, refine_factor, parent_Imin,
 parent_Imax, parent_Jmin, parent_Jmax等）初始化为特定值或空值。

3.如果是3D应用（由is3d判断），那么将一些字段（如Vtransform, Vstretching, theta_s,
 theta_b, Tcline, hc, s_rho, Cs_r, s_w, Cs_w等）初始化为空值。

4.从NetCDF文件中读取spherical字段，并根据其值更新spherical变量。将spherical,uniform,
curvilinear,vector_rotation等字段赋给输出变量Gout。

5.如果Ginp是结构体，代码将检查其是否包含字段spherical。如果包含，那么将其值赋给spherical。

另外，在MATLAB中，||和or都是逻辑运算符，用于表示"或"（OR）的操作，但它们在行为上有一些不同。
||是一个短路逻辑运算符。它会从左到右评估其操作数，一旦找到一个为真（非零），就会立即停止评
估并返回真。这意味着，如果||左边的表达式为真，那么它就不会去评估右边的表达式。

or也是一个逻辑运算符，表示"或"的操作，但它不是短路运算符。无论左边的表达式是否为真，
它都会去评估右边的表达式。

所以，虽然||和or在大多数情况下可以互换使用，但在考虑代码性能以及需要控制代码执行流程
的情况下（例如，当右边的表达式包含函数调用或其他可能产生副作用的操作时），可根据实际
情况选择。
%}
if (~isstruct(Ginp))
  if (any(strcmp(anames,'coarse_factor')) ||                            ...
      any(strcmp(anames,'refine_factor')))
    index = strcmp(anames,'parent_grid');
    if (any(index))
      Gout.parent_grid = Ginfo.Attributes(index).Value;
    end
  end
  Gout.TimeRecord    = TimeRecord;
  Gout.Lm            = [];
  Gout.Mm            = [];
  Gout.N             = [];
  Gout.coarse_factor = 0;
  if (any(strcmp(anames,'coarse_factor')))
    Gout.parent_Imin = [];
    Gout.parent_Imax = [];
    Gout.parent_Jmin = [];
    Gout.parent_Jmax = [];
  end
  Gout.refine_factor = 0;
  Gout.coarse_factor = [];
  if (any(strcmp(anames,'refine_factor')))
    Gout.parent_Imin = [];
    Gout.parent_Imax = [];
    Gout.parent_Jmin = [];
    Gout.parent_Jmax = [];
  end

  if (is3d)
    Gout.Vtransform    = [];
    Gout.Vstretching   = [];
    Gout.theta_s       = [];
    Gout.theta_b       = [];
    Gout.Tcline        = [];
    Gout.hc            = [];
    Gout.s_rho         = [];
    Gout.Cs_r          = [];
    Gout.s_w           = [];
    Gout.Cs_w          = [];
  end
    
  spherical = nc_read(Ginp,'spherical');
  
  if (ischar(spherical))
    if (spherical == 'T' || spherical == 't')
      spherical = 1;
    else
      spherical = 0;
    end
  end
  Gout.spherical = spherical;
  Gout.uniform = false;
  Gout.curvilinear = false;
  Gout.vector_rotation = false;
else
  if (~isempty(Ginp.spherical))
    spherical = Ginp.spherical;
  end
end

%--------------------------------------------------------------------------
% 三.设置变量列表为进程和开关
%--------------------------------------------------------------------------

% 定义前处理变量列表为varhor(horizonal variables)
varhor  = {'x_rho', 'y_rho', 'x_psi', 'y_psi',                          ...
           'x_u', 'y_u',  'x_v', 'y_v'};
%{
以下if语句的逻辑是：根据输入的ROMS网格是否为曲面，往变量列表中加入不同的变量。
%}
if (~spherical)
  varhor = [varhor, 'x_perimeter', 'y_perimeter',                       ...
                    'x_rho_west',  'y_rho_west',                        ...
                    'x_rho_east',  'y_rho_east',                        ...
                    'x_rho_south', 'y_rho_south',                       ...
                    'x_rho_north', 'y_rho_north',                       ...
                    'x_u_west',    'y_u_west',                          ...
                    'x_u_east',    'y_u_east',                          ...
                    'x_u_south',   'y_u_south',                         ...
                    'x_u_north',   'y_u_north',                         ...
                    'x_v_west',    'y_v_west',                          ...
                    'x_v_east',    'y_v_east',                          ...
                    'x_v_south',   'y_v_south',                         ...
                    'x_v_north',   'y_v_north'];
else
  varhor = [varhor, 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',         ...
                    'lon_u', 'lat_u', 'lon_v', 'lat_v',                 ...
                    'lon_perimeter', 'lat_perimeter',                   ...
	            'lon_rho_west',  'lat_rho_west',                    ...
                    'lon_rho_east',  'lat_rho_east',                    ...
                    'lon_rho_south', 'lat_rho_south',                   ...
                    'lon_rho_north', 'lat_rho_north',                   ...
                    'lon_u_west',    'lat_u_west',                      ...
                    'lon_u_east',    'lat_u_east',                      ...
                    'lon_u_south',   'lat_u_south',                     ...
                    'lon_u_north',   'lat_u_north',                     ...
                    'lon_v_west',    'lat_v_west',                      ...
                    'lon_v_east',    'lat_v_east',                      ...
                    'lon_v_south',   'lat_v_south',                     ...
                    'lon_v_north',   'lat_v_north'];
end

% 向前处理变量列表varhor中添加一些新变量，同时构建记录垂直参数的varver(vertical variables)元组。
varhor  = [varhor, 'mask_rho', 'mask_psi', 'mask_u', 'mask_v',          ...
                   'angle', 'pm', 'pn', 'dndx', 'dmde', 'h', 'f'];
varver  = {'Vtransform', 'Vstretching',                                 ...
           'theta_s', 'theta_b', 'Tcline', 'hc',                        ...
           's_rho', 'Cs_r', 's_w', 'Cs_w',                              ...
           'Hz', 'z_rho', 'z_u', 'z_v'};

%{
如果是3d的例子，就将varhor和varver合并成我们需要的变量列表；若不是3d的例子，则我们
的变量列表就不需要包含varver。
%}
if (is3d)
  varlist = [varhor, varver];
else
  varlist = varhor;
end


if (isstruct(Ginp))

  Gout = Ginp;
  
  got.N    = false;
  got.zeta = false;

  for var = varlist
    field = char(var);
    got.(field) = false;
    if (isfield(Gout,field))
      got.(field) = true;
    end
  end

else

  process.horizontal = true;

  got.N              = false;
  got.zeta           = false;
  got.lon_coast      = any(strcmp(vnames,'lon_coast'));
  got.lat_coast      = any(strcmp(vnames,'lat_coast'));

  for var = varlist
    field = char(var);
    got.(field) = false;
  end
   
  for n = 1:NvarsG
    field = char(Ginfo.Variables(n).Name);
    if (isfield(got,field))
      got.(field) = true;
    end

% If "Ginp" is a ROMS output NetCDF file, read in vertical coordinate
% parameters.

    switch field  
      case 'Vtransform'
        Gout.(field) = nc_read(Ginp,field);
      case 'Vstretching'
        Gout.(field) = nc_read(Ginp,field);
      case 'theta_s'
        Gout.(field) = nc_read(Ginp,field);
      case 'theta_b'
        Gout.(field) = nc_read(Ginp,field);
      case 'Tcline'
        Gout.(field) = nc_read(Ginp,field);
      case 'hc'
        Gout.(field) = nc_read(Ginp,field);
      case {'s_rho', 'sc_r'}
        Gout.(field) = nc_read(Ginp,field);
        Gout.N = length(Gout.(field));
        got.N  = true;
      case 'Cs_r'
        Gout.(field) = nc_read(Ginp,field);
      case {'s_w', 'sc_w'}
        Gout.(field) = nc_read(Ginp,field);
      case 'Cs_w'
        Gout.(field) = nc_read(Ginp,field);
      case 'zeta'
        if (is3d)
          process.vertical = true;
        end
        got.zeta  = true;
        zeta_file = Ginp;
    end
  end

% If extracted from finer grid, get coarseness factor.

  got.coarse_factor = any(strcmp(anames,'coarse_factor'));
  if (got.coarse_factor)
    index = strcmp(anames,'coarse_factor');
    Gout.coarse_factor = double(Ginfo.Attributes(index).Value);

    index = strcmp(anames,'parent_Imin');
    Gout.parent_Imin   = double(Ginfo.Attributes(index).Value);
    index = strcmp(anames,'parent_Imax');
    Gout.parent_Imax   = double(Ginfo.Attributes(index).Value);

    index = strcmp(anames,'parent_Jmin');
    Gout.parent_Jmin   = double(Ginfo.Attributes(index).Value);
    index = strcmp(anames,'parent_Jmax');
    Gout.parent_Jmax   = double(Ginfo.Attributes(index).Value);
  end

% If refinement grid, get refinement factor.

  got.refine_factor = any(strcmp(anames,'refine_factor'));
  if (got.refine_factor)
    index = strcmp(anames,'refine_factor');
    Gout.refine_factor = double(Ginfo.Attributes(index).Value);

    index = strcmp(anames,'parent_Imin');
    Gout.parent_Imin   = double(Ginfo.Attributes(index).Value);
    index = strcmp(anames,'parent_Imax');
    Gout.parent_Imax   = double(Ginfo.Attributes(index).Value);

    index = strcmp(anames,'parent_Jmin');
    Gout.parent_Jmin   = double(Ginfo.Attributes(index).Value);
    index = strcmp(anames,'parent_Jmax');
    Gout.parent_Jmax   = double(Ginfo.Attributes(index).Value);
  end

% If vertical parameters "Vtransform" and "Vstretching" are not found
% in primary file "Ginp", set their default values for backward
% compatibility.

  if (process.vertical)  
    if (~got.Vtransform)
      Gout.Vtransform = 1;
      got.Vtransform  = true;   
    end

    if (~got.Vstretching)
      Gout.Vtransform = 1;
      got.Vstretching = true;
    end
  end
  
end

%--------------------------------------------------------------------------
% The "Sinp" argument is a structure array containing application
% vertical coordinate parameters.
%--------------------------------------------------------------------------

if (process.parameters)
  
  parlist = {'N', 'Vtransform', 'Vstretching',                          ...
             'theta_s', 'theta_b', 'Tcline', 'hc'};

  for par = parlist
    field = char(par);
    if (isfield(Sinp,field))
       Gout.(field) = Sinp.(field);
       got.(field)  = true;
    else
      error([' GET_ROMS_GRID: unable to find field "',field,'", in',    ...
             ' in structure: Sinp']);
    end
  end

  process.vertical = true;

end

%--------------------------------------------------------------------------
% The "Sinp" argument is a secondary ROMS output file.  Use the vertical
% coordinates parameters from secodary file to compute depths. Overwrite
% their values in "Gout" structure.
%--------------------------------------------------------------------------

if (process.zeta)
  if (ischar(Sinp))
    Sinfo = nc_inq(Sinp);
    NvarsS = length(Sinfo.Variables);

    got.N           = false;
    got.Vtransform  = false;
    got.Vstretching = false;
    got.theta_s     = false;
    got.theta_b     = false;
    got.Tcline      = false;
    got.hc          = false;
    got.s_rho       = false;
    got.Cs_r        = false;
    got.s_w         = false;
    got.Cs_w        = false;
    got.zeta        = false;

    for n = 1:NvarsS
      field = char(Sinfo.Variables(n).Name);
      switch field
        case 'Vtransform'
          Gout.(field) = nc_read(Sinp,field);
          got.(field)  = true;
        case 'Vstretching'
          Gout.(field) = nc_read(Sinp,field);
          got.(field)  = true;
        case 'theta_s'
          Gout.(field) = nc_read(Sinp,field);
          got.(field)  = true;
        case 'theta_b'
          Gout.(field) = nc_read(Sinp,field);
          got.(field)  = true;
        case 'Tcline'
          Gout.(field) = nc_read(Sinp,field);
          got.(field)  = true;
        case 'hc'
          Gout.(field) = nc_read(Sinp,field);
          got.(field)  = true;
        case {'s_rho', 'sc_r'}
          Gout.(field) = nc_read(Sinp,field);
          got.(field)  = true;
          Gout.N = length(Gout.(field));
          got.N  = true;
        case 'Cs_r'
          Gout.(field) = nc_read(Sinp,field);
          got.(field)  = true;
        case {'s_w', 'sc_w'}
          Gout.(field) = nc_read(Sinp,field);
          got.(field)  = true;
        case 'Cs_w'
          Gout.(field) = nc_read(Sinp,field);
          got.(field)  = true;
        case 'zeta'
          got.(field) = true;
          zeta_file = Sinp;
          process.vertical = true;
      end
    end

% If vertical parameters "Vtransform" and "Vstretching" are not found
% in secondary file, set their default values for backward compatibility.
    
    if (~got.Vtransform)
      Gout.Vtransform  = 1;
      got.Vtransform   = true;
    end

    if (~got.Vstretching)
      Gout.Vstretching = 1;
      got.Vstretching  = true;
    end
    
  end
end

%--------------------------------------------------------------------------
% Add grid variables to structure, except vertical depths.
%--------------------------------------------------------------------------

% If appropriate, read in fields from primary file "Ginp".

if (process.horizontal)
  for var = varhor
    field = char(var);
    if (got.(field))
      Gout.(field) = nc_read(Ginp,field);
    else
      Gout.(field) = [];
    end
  end

  if (got.h)
    [Lr,Mr]=size(Gout.h);
    L = Lr-1;  Lm=L-1;
    M = Mr-1;  Mm=M-1;

    Gout.Lm = Lm;                %  This are horizontal (Lm,Mm) dimensions
    Gout.Mm = Mm;                %  that are specified in ROMS input script
                                 %  RHO-points: Lp = Lm+1,  Mp = Mm+2  

    if (~got.mask_rho)
      Gout.mask_rho = ones(Lr,Mr);
      got.mask_rho  = true;
    end
    if (~got.mask_psi)
      Gout.mask_psi = ones(L,M);
      got.mask_psi  = true;
    end
    if (~got.mask_u)
      Gout.mask_u   = ones(L,Mr);
      got.mask_u    = true;
    end
    if (~got.mask_v)
      Gout.mask_v   = ones(Lr,M);
      got.mask_v    = true;
    end
  
    if (~got.angle)
      Gout.angle = zeros(Lr,Mr);
      got.angle  = true;
    end
  end

% Extract boundary conditions locations.

  if (spherical)
    if (got.lon_rho)
      Gout.lon_rho_west  = Gout.lon_rho(1,:);
      Gout.lon_rho_east  = Gout.lon_rho(end,:);
      Gout.lon_rho_south = Gout.lon_rho(:,1);
      Gout.lon_rho_north = Gout.lon_rho(:,end);
    end
    if (got.lat_rho)
      Gout.lat_rho_west  = Gout.lat_rho(1,:);
      Gout.lat_rho_east  = Gout.lat_rho(end,:);
      Gout.lat_rho_south = Gout.lat_rho(:,1);
      Gout.lat_rho_north = Gout.lat_rho(:,end);
    end

    if (got.lon_u)
      Gout.lon_u_west    = Gout.lon_u(1,:);
      Gout.lon_u_east    = Gout.lon_u(end,:);
      Gout.lon_u_south   = Gout.lon_u(:,1);
      Gout.lon_u_north   = Gout.lon_u(:,end);
    end
    if (got.lat_u)
      Gout.lat_u_west    = Gout.lat_u(1,:);
      Gout.lat_u_east    = Gout.lat_u(end,:);
      Gout.lat_u_south   = Gout.lat_u(:,1);
      Gout.lat_u_north   = Gout.lat_u(:,end);
    end

    if (got.lon_v)
      Gout.lon_v_west    = Gout.lon_v(1,:);
      Gout.lon_v_east    = Gout.lon_v(end,:);
      Gout.lon_v_south   = Gout.lon_v(:,1);
      Gout.lon_v_north   = Gout.lon_v(:,end);
    end
    if (got.lat_v)
      Gout.lat_v_west    = Gout.lat_v(1,:);
      Gout.lat_v_east    = Gout.lat_v(end,:);
      Gout.lat_v_south   = Gout.lat_v(:,1);
      Gout.lat_v_north   = Gout.lat_v(:,end);
    end
  else
    if (got.x_rho)
      Gout.x_rho_west  = Gout.x_rho(1,:);
      Gout.x_rho_east  = Gout.x_rho(end,:);
      Gout.x_rho_south = Gout.x_rho(:,1);
      Gout.x_rho_north = Gout.x_rho(:,end);
    end
    if (got.y_rho)
      Gout.y_rho_west  = Gout.y_rho(1,:);
      Gout.y_rho_east  = Gout.y_rho(end,:);
      Gout.y_rho_south = Gout.y_rho(:,1);
      Gout.y_rho_north = Gout.y_rho(:,end);
    end

    if (got.x_u)
      Gout.x_u_west    = Gout.x_u(1,:);
      Gout.x_u_east    = Gout.x_u(end,:);
      Gout.x_u_south   = Gout.x_u(:,1);
      Gout.x_u_north   = Gout.x_u(:,end);
    end
    if (got.y_u)
      Gout.y_u_west    = Gout.y_u(1,:);
      Gout.y_u_east    = Gout.y_u(end,:);
      Gout.y_u_south   = Gout.y_u(:,1);
      Gout.y_u_north   = Gout.y_u(:,end);
    end

    if (got.x_v)
      Gout.x_v_west    = Gout.x_v(1,:);
      Gout.x_v_east    = Gout.x_v(end,:);
      Gout.x_v_south   = Gout.x_v(:,1);
      Gout.x_v_north   = Gout.x_v(:,end);
    end
    if (got.y_v)
      Gout.y_v_west    = Gout.y_v(1,:);
      Gout.y_v_east    = Gout.y_v(end,:);
      Gout.y_v_south   = Gout.y_v(:,1);
      Gout.y_v_north   = Gout.y_v(:,end);
    end
  end

% Extract perimeter.

  if (spherical)
    if (got.lon_psi && got.lat_psi)
      [L,M]=size(Gout.lon_psi);
      Gout.lon_perimeter = [squeeze(Gout.lon_psi(1:L,1));               ...
                            squeeze(Gout.lon_psi(L,2:M))';              ...
                            squeeze(flipud(Gout.lon_psi(1:Lm,M)));      ...
                            squeeze(fliplr(Gout.lon_psi(1,1:Mm)))'];

      Gout.lat_perimeter = [squeeze(Gout.lat_psi(1:L,1));               ...
                            squeeze(Gout.lat_psi(L,2:M))';              ...
                            squeeze(flipud(Gout.lat_psi(1:Lm,M)));      ...
                            squeeze(fliplr(Gout.lat_psi(1,1:Mm)))'];
    end
  else
    if (got.x_psi && got.x_psi)
      [L,M]=size(Gout.x_psi);
      Gout.x_perimeter = [squeeze(Gout.x_psi(1:L,1));                   ...
                          squeeze(Gout.x_psi(L,2:M))';                  ...
                          squeeze(flipud(Gout.x_psi(1:Lm,M)));          ...
                          squeeze(fliplr(Gout.x_psi(1,1:Mm)))'];

      Gout.y_perimeter = [squeeze(Gout.y_psi(1:L,1));                   ...
                          squeeze(Gout.y_psi(L,2:M))';                  ...
                          squeeze(flipud(Gout.y_psi(1:Lm,M)));          ...
                          squeeze(fliplr(Gout.y_psi(1,1:Mm)))'];
    end
  end
  
% Determine "uniform", "curvilinear", and "vector_rotation" switches.

  if (got.pm && got.pn)
    if (length(unique(Gout.pm(:))) == 1 &&                              ...
        length(unique(Gout.pn(:))) == 1)
      Gout.uniform = true;
    end

    if (length(unique(Gout.pm(:))) > 1 ||                               ...
        length(unique(Gout.pn(:))) > 1)
      Gout.curvilinear = true;
    end
  end
  
  if (got.angle)
    if (length(unique(Gout.angle(:))) > 1 ||                            ...
               unique(Gout.angle(:))  > 0)
      Gout.vector_rotation = true;
    end
  end

% If available, process coastline data.

  if (got.lon_coast)
    Gout.lon_coast = nc_read(Ginp,'lon_coast');
    Gout.lat_coast = nc_read(Ginp,'lat_coast');
  end

end

%--------------------------------------------------------------------------
% Add vertical depths.
%--------------------------------------------------------------------------

if (process.vertical && is3d)

  if (~got.Vtransform  || isempty(Gout.Vtransform ))
    error([' GET_ROMS_GRID: unassigned field ''Vtransform''',           ...
           ' in structure: Gout']);
  end
  if (~got.Vstretching || isempty(Gout.Vstretching))
    error([' GET_ROMS_GRID: unassigned field ''Vstretching''',          ...
           ' in structure: Gout']);
  end
  if (~got.theta_s     || isempty(Gout.theta_s    ))
    error([' GET_ROMS_GRID: unassigned field ''theta_s''',              ...
           ' in structure: Gout']);
  end
  if (~got.theta_b     || isempty(Gout.theta_b    ))
    error([' GET_ROMS_GRID: unassigned field ''theta_b''',              ...
           ' in structure: Gout']);
  end
  if (~got.hc          || isempty(Gout.hc         ))
    error([' GET_ROMS_GRID: unassigned field ''hc''',                   ...
           ' in structure: Gout']);
  end
  if (~got.N           || isempty(Gout.N          ))
    error([' GET_ROMS_GRID: unassigned field ''N''',                    ...
           ' in structure: Gout']);
  end      

  h = Gout.h;
  if (~isempty(TimeRecord))
    zeta = nc_read(zeta_file,'zeta',TimeRecord);
  else
    zeta = zeros(size(h));
    Gout.TimeRecord = 'Computing unperturbed depths, zeta=0';
  end
   
  if (isempty(h))
    disp(' ')
    disp('   GET_ROMS_GRID - input file does not have grid data:');
    disp(['            Ginp:  ',Ginp]);
  else
    if (~(got.s_rho || got.Cs_w))
      kgrid = 0;
      [Gout.s_rho,Gout.Cs_r] = stretching(Gout.Vstretching,             ...
                                          Gout.theta_s,                 ...
                                          Gout.theta_b,                 ...
                                          Gout.hc, Gout.N,              ...
                                          kgrid, false);
    end

    if (~(got.s_w || got.Cs_w))
      kgrid = 1;
      [Gout.s_w,  Gout.Cs_w] = stretching(Gout.Vstretching,             ...
                                          Gout.theta_s,                 ...
                                          Gout.theta_b,                 ...
                                          Gout.hc, Gout.N, kgrid);
    end

    igrid = 1;
    Gout.z_r = set_depth(Gout.Vtransform, Gout.Vstretching,             ...
                         Gout.theta_s, Gout.theta_b, Gout.hc,           ...
                         Gout.N, igrid, h, zeta, false);

    igrid = 3;
    Gout.z_u = set_depth(Gout.Vtransform, Gout.Vstretching,             ...
                         Gout.theta_s, Gout.theta_b, Gout.hc,           ...
                         Gout.N, igrid, h, zeta, false);

    igrid = 4;
    Gout.z_v = set_depth(Gout.Vtransform, Gout.Vstretching,             ...
                         Gout.theta_s, Gout.theta_b, Gout.hc,           ...
                         Gout.N, igrid, h, zeta, false);

    igrid = 5;
    Gout.z_w = set_depth(Gout.Vtransform, Gout.Vstretching,             ...
                         Gout.theta_s, Gout.theta_b, Gout.hc,           ...
                         Gout.N, igrid, h, zeta, false);

    N = Gout.N;
    Gout.Hz = Gout.z_w(:,:,2:N+1) - Gout.z_w(:,:,1:N);
  end

end

return