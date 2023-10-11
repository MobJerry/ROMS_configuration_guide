% jcw, revised, Feb 10, 2019
% modified by pzy,Sep 27,2023

%{
get_ijrg函数是用于找到我们的ROMS自定义网格在hycom网格中的索引，从而确
定hycom网格的一个子集来获取数据。

该函数的输入参数有：1.数据下载网址的url；2.网格路径名(modelgrid)；3.theta_s, 
 theta_b,Tcline, N, Vtransform, Vstretching为输入网格的一些参数。

该函数的输出有：1.gn为包含给定网格所有变量的结构数组；2.clm也是一个结构体，里面
包含了roms_master_climatology_coawst_mw生成clm文件所需的hycom网格索引信息、插值
深度等。
%}

function [gn, clm]=get_ijrg(url, modelgrid, theta_s, theta_b, Tcline, N, Vtransform, Vstretching)

% 首先，读取ROMS网格信息

disp('正在获取ROMS网格各维度数据...');

% 下面，定义一个结构体变量sinp，用于存放roms网格的信息

Sinp.theta_s     = theta_s;      % theta_s是表层控制参数
Sinp.theta_b     = theta_b;      % theta_b是底层控制参数
Sinp.Tcline      = Tcline;       % Tcline是表/底拉伸宽度
Sinp.N           = N;            % N是垂直层级的数量
Sinp.Vtransform  = Vtransform;   % Vtransform是网格对应的垂直转换算法
Sinp.Vstretching = Vstretching;  % Vstretching是网格对应的垂直拉伸算法

if (Vtransform==1)               % 该段if语句主要是为了确定hc的值。详见wikiROMS
  h = ncread(modelgrid,'h');
  hmin=min(h(:));
  hc=min(max(hmin,0),Tcline);
elseif (Vtransform==2)
  hc=Tcline;
end

Sinp.hc          = hc;               % hc为ROMS中所使用的临界深度
gn = get_roms_grid(modelgrid,Sinp);  % 使用get_roms_grid函数获取包含自定义网格
                                     % 所有变量的结构数组gn

%{
这里用shiftdim函数，将gn中的某些变量维度向左移动两个单位。这些变量的维度
为L*M*N，移动后变为N*L*M。
%}
gn.z_r = shiftdim(gn.z_r,2);         
gn.z_u = shiftdim(gn.z_u,2);
gn.z_v = shiftdim(gn.z_v,2);
gn.z_w = shiftdim(gn.z_w,2);

% 读出HYCOM的lon、lat和depth变量

disp(['从',url,'处获取HYCOM网格数据...']);

%{
接下来的四行代码已打上注释。HYCOM的最新数据集(2018-12至今)已没有X、Y变量
，因此建议在下载数据之前先打开对应网址看看变量名称。
%}
% numX = ncread(url,'X');              
% numY = ncread(url,'Y');              
% hycom_lon = ncread(url,'Longitude',[1 1],[length(numX) 1]);% 
% hycom_lat=ncread(url,'Latitude',[1 1],[1 length(numY)]);
hycom_lon = ncread(url,'lon');
hycom_lat = ncread(url,'lat');
hycom_depth = ncread(url,'depth');

% 获取ROMS网格的经纬度范围(xl:x_left,xr:x_right,yb:y_bottom,yt:y_top)

disp('获取ROMS网格维度...');
xl = min(min(gn.lon_rho)); xr = max(max(gn.lon_rho));
yb = min(min(gn.lat_rho)); yt = max(max(gn.lat_rho));

% 接下来调整ROMS的网格范围落在HYCOM网格之中

disp('优化网格维度...');
xg = hycom_lon;
xg(xg >= 180) = (xg(xg >= 180) - 360);   % 将原本范围0-360°调整为-180-180°
yg = hycom_lat';                          % 纬度范围保持-80-90°不变

% 利用find函数将包含在hycom网格(xg yg)之中的roms网格(xl xr yb yt)索引找出来

[ii] = find(xg >= xl & xg <= xr);
[jj] = find(yg >= yb & yg <= yt);

% 在ROMS网格的边界处多取一个值

ig0=(min(ii)-1); ig1=(max(ii)+1); jg0=(min(jj)-1); jg1=(max(jj)+1);

% 做一些严谨的处理，获得完全落在hycom网格中的ROMS网格的索引

ig0 = max(ig0, 1);
jg0 = max(jg0, 1);
ig1 = min(ig1, length(hycom_lat));
jg1 = min(jg1, length(hycom_lon));

% 把这些索引保存为字符串

irg2=[num2str(ig0) ':' num2str(ig1)];
jrg2=[num2str(jg0) ':' num2str(jg1)];

% 最后生成对应的clm文件

clm.lon=double(xg(ig0:ig1));
clm.lat=double(yg(jg0:jg1));
clm.z=double(hycom_depth);
clm.irg2=irg2;
clm.jrg2=jrg2;
clm.ig0=ig0;
clm.ig1=ig1;
clm.jg0=jg0;
clm.jg1=jg1;

