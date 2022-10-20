# SegdProcessMasw
# 脚本方法主要参考Park et al. 1998.，其中source文件夹为segy数据，运行脚本main.py
######
# segd_analysis.get_dispersion(st_new,dx,cmin,cmax,dc,fmax)
# st_new 道集数据，以obspy读取数据为准
# dx 道间距
# cmin,cmax 最小，最大相速度m/s
# dc 相速度采样间隔m/s
# fmax 最大频率
#######
# segd_analysis.source_xkm(ShotGather,SourceTrace,CutTraceNum)
# 截取震源附近的道集仅适用于segd文件
# ShotGather 道集数据，以obspy读取的segd文件为准，需要包含完整的检波点坐标与震源坐标数据
# SourceTrace 新截取的道集中的震源位置
# CutTraceNum 新截取的道集中总道数
#######




###
reference
The phase shift method described in: Park, C.B., Miller, R.D. and Xia, J., 1998, January. Imaging dispersion curves of surface waves on multi-channel record. In 1998 SEG Annual Meeting. Society of Exploration Geophysicists.
https://github.com/luan-th-nguyen/PyDispersion
