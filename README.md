# ROMS_configuration_guide
该项目包含了对用户自定义算例进行ROMS一些必要的输入文件的制作教程，以及头文件选项的选择方法。持续更新。

这里，mfile是一个matlab工具包，用于生成ROMS网格所需的netcdf input文件。只需将mfile添加到自己的matlab
路径中，然后打开roms_clm/roms_master_climatology_coawst_mw![image](https://github.com/MobJerry/ROMS_configuration_guide/assets/94905594/6e02edf3-5a22-40f1-977d-a257e6e7df00)
完成用户输入，运行便可得到ini、bdy和clm文件。

注意，代码需要根据选用数据集的不同而要进行适当的更改。我已为代码添加详细的中文注释，按照注释更改即可。
此外，为了性能上的提升和便于用户阅读和维护，以后会不时更改代码。

