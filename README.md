# ROMS_configuration_guide

如果你也对网上ROMS教程太少而苦恼，如果你也对ROMSwiki上的一堆英文看不懂而发愁，如果你也对ROMS不能迟迟
入门而沮丧，来这里看看吧！

该项目包含了对用户自定义算例进行ROMS一些必要的输入文件的制作教程，以及头文件选项的选择方法。持续更新。

## 目录

- [Security](#security)

该目录下包含ROMS网格生成应用GridBuilder的几个实用教程。

这里，mfile是一个matlab工具包，用于生成ROMS网格所需的netcdf input文件。只需将mfile添加到自己的matlab
路径中，然后打开roms_clm/roms_master_climatology_coawst_mw完成用户输入，运行便可得到ini、bdy和clm文件。
mfile中有一个文件夹rutgers太大，无法直接上传；于是我将其压缩后上传，大家下载mfile到本地之后记得解压缩
rutgers就可以（刚玩github不熟操作，菜鸟一个，大家轻喷！）

注意，代码需要根据选用数据集的不同而要进行适当的更改。我已为代码添加详细的中文注释，按照注释更改即可。
此外，为了性能上的提升和便于用户阅读和维护，以后会不时更改代码。源码并非我所写，我只是做出适当更改，以及
添加中文注释方便大家查看，因此仅可做学习交流用途，切不可用于盈利。如若发现，依法追究法律责任。

