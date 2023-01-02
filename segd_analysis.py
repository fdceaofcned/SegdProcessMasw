import obspy
from numpy import matrix as mat
import numpy as np
from matplotlib import pyplot as plt
from read_segd import read_segd # 读取segd文件
from dispersion import get_dispersion
import time, os, sys
# from pysurf96 import surf96 # 正演频散曲线此处需要conda install libpython库 编译要mgwin64库手动下载安装
# 地震数据操作
def connon_psition(st0): # 提取震源位置
    cannon_e = st0[0].stats.segd.source_easting
    cannon_n = st0[0].stats.segd.source_northing
    cannon_h = st0[0].stats.segd.source_elevation
    return cannon_e,cannon_n,cannon_h
def plot_segd_plus(st0): # 绘图终极版
    x,y,z,e,n,h = segd_position(st0)
    plot_segd(x, y, e, n)
    return plt.show()
def segd_position(st0): # 提取检波点坐标
    station_position_e, station_position_n, station_position_h = [], [], []
    for i in range(0,st0.count(),1):
        station_position_e.append(st0[i].stats.segd.receiver_point_easting) # 东西坐标
        station_position_n.append(st0[i].stats.segd.receiver_point_northing) # 南北坐标
        station_position_h.append(st0[i].stats.segd.receiver_point_elevation) # 高程
    cannon_e,cannon_n,cannon_h = connon_psition(st0)
    return np.array(station_position_e), np.array(station_position_n), np.array(station_position_h), cannon_e, cannon_n, cannon_h
def rec_position(st0): # 提取检波点坐标
    station_position_e, station_position_n, station_position_h = [], [], []
    for i in range(0,st0.count(),1):
        station_position_e.append(st0[i].stats.segd.receiver_point_easting) # 东西坐标
        station_position_n.append(st0[i].stats.segd.receiver_point_northing) # 南北坐标
        station_position_h.append(st0[i].stats.segd.receiver_point_elevation) # 高程
    return station_position_e,station_position_n,station_position_h
# 绘制平面点位分布带序号标定
def plot_segd(position_e, position_n, cannon_e, cannon_n):
    del_e = np.where(position_e==0)
    position_e = np.append(np.delete(position_e,del_e),cannon_e)
    position_n = np.append(np.delete(position_n,del_e),cannon_n)
    colors = np.zeros(len(position_e)-1)+32
    colors = np.append(colors,64)
    n = np.arange(len(position_e)) + 1
    fig,ax = plt.subplots()
    ax.scatter(position_e,position_n,c=colors)
    for i,txt in enumerate(n):
        ax.annotate(txt,(position_e[i],position_n[i]))
    return plt.show()
# 绘制平面点位分布
def plot_segd_old(position_e, position_n, cannon_e, cannon_n):
    del_e = np.where(position_e==0)
    position_e = np.append(np.delete(position_e,del_e),cannon_e)
    position_n = np.append(np.delete(position_n,del_e),cannon_n)
    colors = np.zeros(len(position_e)-1)+32
    colors = np.append(colors,64)
    plt.scatter(position_n,position_e,c=colors)
    return plt.show()
# 提取数据矩阵 输入起点，终点，文件，采样率
def rayleigh_data(start_point,end_point,st0,rate): 
    data_array,time_array = np.array([]),np.array([])
    for i in range(start_point,end_point+1,1):
        tr = st0[i]
        tr.decimate(rate)
        print(len(tr.data))
        data_array = np.append(data_array,tr.data)
    return np.reshape(data_array,(end_point-start_point+1,len(tr.data)))
import linecache
import numpy as np
def get_line(file, nums_line): # 封装读取txt行函数（文件，指定行）
    line_aim = linecache.getline(file, nums_line).strip()
    return line_aim
# 频散能量图绘制
def draw_fk(im,ax,f,c,img,fmax_idx):
    timestamp = time.strftime(".%Y.%m.%d-%H.%M.%S", time.localtime())
    ax.imshow(img[:,:],aspect='auto',origin='lower',extent=(f[0],f[fmax_idx],c[0],c[-1]),interpolation='bicubic',cmap='jet')
    ax.grid(linestyle='--',linewidth=0.5)
    ax.set_xlabel('Frequency [Hz]', fontsize=10)
    ax.set_ylabel('Phase velocity [m/s]', fontsize=10)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 10)
    return im,timestamp
# 频散能量图等值线效果待优化
def draw_fc_contour(plus_img,fmax,cmin,cmax):
    import numpy as np
    f_index = np.linspace(0,fmax,len(plus_img))
    c_index = np.linspace(cmin,cmax,len(plus_img[0]))
    # f_index = np.arange(0,len(plus_img)) # c_index = np.arange(0,len(plus_img[0]))
    f_mesh,c_mesh = np.meshgrid(f_index,c_index)
    C = plt.contour(c_index,f_index,plus_img,18,colors='black',linewidth=.175)
    plt.contourf(c_index,f_index,plus_img,18)
    return C
# 绘制截取图像
def transformer_plus(source,target,x_lat,y_lon): # 批量list转换坐标,(源wkid，目标wkid，x或纬度，y或经度）
    from pyproj import Transformer
    target_x,target_y = [],[]
    transformer = Transformer.from_crs(source,target)
    for i in range(0,len(x_lat),1):
        per_x,per_y = transformer.transform(x_lat[i],y_lon[i])
        target_x.append(per_x)
        target_y.append(per_y)
    return target_x,target_y
def source_xkm(st0,fix,offset): #将炮点附近点提取出来（炮集，新道集相对震源位置，偏移道数）
    e_point_list, n_point_list, h_point_list, e_source, n_source, h_source = segd_position(st0)
    x_test = (e_source - e_point_list)**2
    y_test = (n_source - n_point_list)**2
    d_en = x_test + y_test # 遍历震源与各点距离
    min_point = np.argmin(d_en) # 找出离震源最近的点位
    if min_point > fix:
        start_point = min_point - fix
    else:
        start_point = 0
    # divide_trace(1,num,start_point,st0,path,va_save_path,file_name)
    st_new = st0[start_point:start_point+offset:1]
    return st_new
def calc_base_parament(st_new): # (traces）
    import numpy as np
    e_point_list, n_point_list, h_point_list, e_source, n_source, h_source = segd_position(st_new)
    z1 = np.polyfit(e_point_list, n_point_list, 1)
    k = z1[0]
    b = z1[1]
    k2 = -(1/k)
    b2 = -k2 * e_source + n_source
    x_new = (b - b2)/(k2 - k)
    y_new = (k2*b - k*b2)/(k2 - k)
    # p_min = np.array([e_point_list[0],n_point_list[0]])
    # p_max = np.array([e_point_list[-1],n_point_list[-1]])
    cross_point = np.array([x_new,y_new])
    # d_first = np.linalg.norm(p_min - cross_point)
    # d_last = np.linalg.norm(p_max - cross_point)
    # d_line = np.linalg.norm(p_min - p_max)
    # dx = d_line/(len(st_new) - 1)
    return cross_point
def plot_sesmic_wave(st_new,zone,record_len,count_tr,imgsize): # for segd type file only imgsize = (7.3,4.1)
    import utm
    from obspy.geodetics import gps2dist_azimuth
    from matplotlib.transforms import blended_transform_factory
    cross_point = calc_base_parament(st_new)
    eq_lat = utm.to_latlon(cross_point[0],cross_point[1],zone,'R')[0]
    eq_lon = utm.to_latlon(cross_point[0],cross_point[1],zone,'R')[1]
    i = 0
    for tr in st_new:
        i = i + 1
        receive_lat = tr.stats.segd.receiver_point_northing
        receive_lon = tr.stats.segd.receiver_point_easting
        point_lalo = utm.to_latlon(receive_lon,receive_lat,zone,'R')
        if eq_lon - point_lalo[1] > 0:
            pos = -1
        else:
            pos = 1
        tr.stats.distance = pos*gps2dist_azimuth(point_lalo[0],point_lalo[1],eq_lat,eq_lon)[0]
        tr.stats.network = str(i)
    fig = plt.figure(figsize = imgsize)
    st_new.plot(type='section',plot_dx=1e3,recordlength=record_len,time_down=True,linewidth=.5,grid_linewidth=.25,show=False,fig=fig)
    ax = fig.axes[0]
    transform = blended_transform_factory(ax.transData, ax.transAxes)
    for i in range(0,len(st_new),count_tr):
        tr = st_new[i]
        ax.text(tr.stats.distance/1e3,1.0,int(tr.stats.network)-1,rotation=0,va="bottom",ha="center",transform=transform,zorder=10)
    return fig
# 坐标分区标签
def cut_str(str_aim,key): #
    res = list(filter(None,str_aim.split(key))) #
    return res
### image plot ###
def gauss_func(x,mu,sigma,a):
    # y = (((2*math.pi)**(-0.5))*(sigma**-1))*((math.e**-((x-mu)**2)/(2*sigma**2)))**0.5
    y = a*(np.exp(-((x - mu)**2)/(2*sigma**2)) / (sigma * np.sqrt(2*np.pi)))
    return y
def SmoothLine(x,y,multiple): # scatter line smooth. multiple is the number for smooth level
    from scipy.interpolate import interp1d
    new_inter = interp1d(x,y,kind='cubic')
    new_x = np.linspace(x[0],x[-1],multiple*len(x))
    new_y = new_inter(new_x)
    return new_x,new_y
def LimitFreq(x,y,limit):
    lim_x = [i for i in x if i < limit]
    end = len(lim_x)
    lim_y = y[0:end:1]
    return lim_x,lim_y
def AddDcLine(file_path,ax,multiple,limit): # multiple is the number for smooth level
    check_list = np.loadtxt(file_path)
    freq = check_list[:,0]
    vr = check_list[:,1]
    # freq = list(map(float,freq))
    # vs = list(map(float,vs))
    freq,vr = SmoothLine(freq,vr,multiple)
    freq,vr = LimitFreq(freq,vr,limit)
    ax.plot(freq,vr)
    return ax
