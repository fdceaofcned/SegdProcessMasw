import obspy
from matplotlib import pyplot as plt
import segd_analysis
### one shot-gather data MASW ### 
ShotGatherPath = 'source/00000001.00000192.segd.sgy'
st_new = obspy.read(ShotGatherPath)
f,c,img,fmax_idx,U,t = segd_analysis.get_dispersion(st_new,40,10.,4000.0,20.,12)
img = img / img.max()
img = 1 - segd_analysis.gauss_func(img,0,1,1)
img = img / img.max()
im,ax = plt.subplots(figsize=(8,4.5))
fvr_pic,timestamp = segd_analysis.draw_fk(im,ax,f,c,img,fmax_idx)
fvr_pic.savefig('00000001.00000192.svg')
### end ###

### folder process ###
segy_path = 'source/'
save_img = ''
masw_para = [40,10.,4000.0,20.,12] # [道间距,最小瑞雷波速,最大瑞雷波速,波速采样,最大频率]
gauss_opt = [0,2,1] # [mu,sigma,a]
plotsize = (8,4.5)
segd_analysis.F_vr_folder(save_img,segy_path,masw_para,gauss_opt,plotsize)
### end ###

### one shot-gather data MASW for segd only ###
# ShotGatherPath = 'source_segd/00000001.00000192.segd'
# st = segd_analysis.read_segd(ShotGatherPath)
# st_new = segd_analysis.source_xkm(st,36,72) # cut the nearly trace of source for segd data (ShotGather,SourceTrace,CutTraceNum)
import pickle
save_path = 'source_segd/segd_pack.dump' # 因为segd文件保密，所以只上传了一个内存转储的文件
segdfile = open(save_path,'rb')
st_new = pickle.load(segdfile)
segdfile.close()
f,c,img,fmax_idx,U,t = segd_analysis.get_dispersion(st_new,40,10.,4000.0,20.,12)
img = img / img.max()
img = 1 - segd_analysis.gauss_func(img,0,1,1)
img = img / img.max()
im,ax = plt.subplots(figsize=(8,4.5))
fvr_pic,timestamp = segd_analysis.draw_fk(im,ax,f,c,img,fmax_idx)
fvr_pic.savefig('00000001.00000192.svg')
### end ###

### folder process ###
segy_path = 'source_segd/'
save_img = ''
masw_para = [40,10.,4000.0,20.,12] # [道间距,最小瑞雷波速,最大瑞雷波速,波速采样,最大频率]
gauss_opt = [0,2,1] # [mu,sigma,a]
plotsize = (8,4.5)
cut_source = [36,72] # [SourceTrace,CutTraceNum]
segd_analysis.F_vr_segd(save_img,segy_path,cut_source,masw_para,gauss_opt,plotsize)
### end ###
