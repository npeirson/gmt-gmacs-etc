import os
import json
import pandas as pd
import numpy as np
import config as cfg
import defaults as dfs
import values as edl
from bokeh.models import ColumnDataSource


# jsonify_data = lambda data,filename: json.dump(data,filename) # for an optimization in later versions

''' import data '''
# explicit... not very elegant, but pandas don't like loops
sb1_x,sb1_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[0]),usecols=(0,1),unpack=True)
sb2_x,sb2_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[1]),usecols=(0,1),unpack=True)
sb3_x,sb3_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[2]),usecols=(0,1),unpack=True)
sb4_x,sb4_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[3]),usecols=(0,1),unpack=True)
sb5_x,sb5_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[4]),usecols=(0,1),unpack=True)
sb6_x,sb6_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[5]),usecols=(0,1),unpack=True)
s0_x,s0_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[6]),usecols=(0,1),unpack=True)
sa_x,sa_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[7]),usecols=(0,1),unpack=True)
sb_x,sb_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[8]),usecols=(0,1),unpack=True)
sc_x,sc_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[9]),usecols=(0,1),unpack=True)
bulge_x,bulge_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[10]),usecols=(0,1),unpack=True)
ellipticals_x,ellipticals_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[11]),usecols=(0,1),unpack=True)
lbg_all_flam_x,lbg_all_flam_y = np.loadtxt(os.path.join(edl.galaxy_path,edl.galaxy_files[12]),usecols=(0,1),unpack=True)

o5v_x,o5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[0]),usecols=(0,1),unpack=True)
b0v_x,b0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[1]),usecols=(0,1),unpack=True)
b57v_x,b57v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[2]),usecols=(0,1),unpack=True)
a0v_x,a0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[3]),usecols=(0,1),unpack=True)
a5v_x,a5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[4]),usecols=(0,1),unpack=True)
f0v_x,f0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[5]),usecols=(0,1),unpack=True)
f5v_x,f5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[12]),usecols=(0,1),unpack=True)
g0v_x,g0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[6]),usecols=(0,1),unpack=True)
g5v_x,g5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[7]),usecols=(0,1),unpack=True)
k0v_x,k0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[8]),usecols=(0,1),unpack=True)
k5v_x,k5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[9]),usecols=(0,1),unpack=True)
m0v_x,m0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[10]),usecols=(0,1),unpack=True)
m5v_x,m5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[11]),usecols=(0,1),unpack=True)

filter_photonux_x,filter_photonux_y = np.loadtxt(os.path.join(edl.filter_path,edl.filter_files[0]),usecols=(0,1),unpack=True)
filter_photonb_x,filter_photonb_y = np.loadtxt(os.path.join(edl.filter_path,edl.filter_files[1]),usecols=(0,1),unpack=True)
filter_photonv_x,filter_photonv_y = np.loadtxt(os.path.join(edl.filter_path,edl.filter_files[2]),usecols=(0,1),unpack=True)
filter_photonr_x,filter_photonr_y = np.loadtxt(os.path.join(edl.filter_path,edl.filter_files[3]),usecols=(0,1),unpack=True)
filter_photoni_x,filter_photoni_y = np.loadtxt(os.path.join(edl.filter_path,edl.filter_files[4]),usecols=(0,1),unpack=True)
filter_u_x,filter_u_y = np.loadtxt(os.path.join(edl.filter_path,edl.filter_files[5]),usecols=(0,1),unpack=True,delimiter=',')
filter_g_x,filter_g_y = np.loadtxt(os.path.join(edl.filter_path,edl.filter_files[6]),usecols=(0,1),unpack=True,delimiter=',')
filter_r_x,filter_r_y = np.loadtxt(os.path.join(edl.filter_path,edl.filter_files[7]),usecols=(0,1),unpack=True,delimiter=',')
filter_i_x,filter_i_y = np.loadtxt(os.path.join(edl.filter_path,edl.filter_files[8]),usecols=(0,1),unpack=True,delimiter=',')
filter_z_x,filter_z_y = np.loadtxt(os.path.join(edl.filter_path,edl.filter_files[9]),usecols=(0,1),unpack=True,delimiter=',')

skyfile_00d_x,skyfile_00d_y = np.loadtxt(os.path.join(edl.skyfiles_path,edl.skyfiles[0]),usecols=(0,1),unpack=True)
skyfile_03d_x,skyfile_03d_y = np.loadtxt(os.path.join(edl.skyfiles_path,edl.skyfiles[1]),usecols=(0,1),unpack=True)
skyfile_07d_x,skyfile_07d_y = np.loadtxt(os.path.join(edl.skyfiles_path,edl.skyfiles[2]),usecols=(0,1),unpack=True)
skyfile_10d_x,skyfile_10d_y = np.loadtxt(os.path.join(edl.skyfiles_path,edl.skyfiles[3]),usecols=(0,1),unpack=True)
skyfile_14d_x,skyfile_14d_y = np.loadtxt(os.path.join(edl.skyfiles_path,edl.skyfiles[4]),usecols=(0,1),unpack=True)

grating_red_x,grating_red_y = np.loadtxt(edl.grating_path[0],usecols=(0,1),unpack=True)
grating_blue_x,grating_blue_y = np.loadtxt(edl.grating_path[1],usecols=(0,1),unpack=True)
ccd_red_x,ccd_red_y = np.loadtxt(edl.ccd_path[0],usecols=(0,1),unpack=True)
ccd_blue_x,ccd_blue_y = np.loadtxt(edl.ccd_path[1],usecols=(0,3),unpack=True)

_dichroic_x,dichroic_y1,dichroic_y2 = np.loadtxt(edl.dichroic_path,usecols=(0,1,2),unpack=True)
dichroic_x = _dichroic_x * 10
atmo_ext_x,atmo_ext_y = np.loadtxt(edl.atmo_ext_path,usecols=(0,1),unpack=True,delimiter=',')

mirror_file = np.fliplr(np.loadtxt(edl.mirror_file,usecols=(0,1),unpack=True,delimiter=' '))


''' coalesced data '''

galaxyfiles = [[sb1_x,sb1_y],[sb2_x,sb2_y],[sb3_x,sb3_y],[sb4_x,sb4_y],[sb5_x,sb5_y],
				[sb6_x,sb6_y],[s0_x,s0_y],[sa_x,sa_y],[sb_x,sb_y],[sc_x,sc_y],
				[bulge_x,bulge_y],[ellipticals_x,ellipticals_y],[lbg_all_flam_x,lbg_all_flam_y]]
starfiles = [[o5v_x,o5v_y],[b0v_x,b0v_y],[b57v_x,b57v_y],
			[a0v_x,a0v_y],[a5v_x,a5v_y],[f0v_x,f0v_y],[f5v_x,f5v_y],
			[g0v_x,g0v_y],[g5v_x,g5v_y],[k0v_x,k0v_y],[k5v_x,k5v_y],
			[m0v_x,m0v_y],[m5v_x,m5v_y]]
filterfiles = [[filter_photonux_x,filter_photonux_y],[filter_photonb_x,filter_photonb_y],
				[filter_photonv_x,filter_photonv_y],[filter_photonr_x,filter_photonr_y],
				[filter_photoni_x,filter_photoni_y],[filter_u_x,filter_u_y],[filter_g_x,filter_g_y],
				[filter_r_x,filter_r_y],[filter_i_x,filter_i_y],[filter_z_x,filter_z_y]]
skyfiles = [[skyfile_00d_x,skyfile_00d_y],[skyfile_03d_x,skyfile_03d_y],[skyfile_07d_x,skyfile_07d_y],
			[skyfile_10d_x,skyfile_10d_y],[skyfile_14d_x,skyfile_14d_y]]
dichroic_x,dichroic_y1,dichroic_y2 = dichroic_x,dichroic_y1,dichroic_y2
grating1,grating2 = [grating_blue_x,grating_blue_y],[grating_red_x,grating_red_y]
ccd1,ccd2 = [ccd_blue_x,ccd_blue_y],[ccd_red_x,ccd_red_y]
atmo_ext_x,atmo_ext_y = atmo_ext_x,atmo_ext_y
mirror_file_x,mirror_file_y = mirror_file[0]*10,mirror_file[1]


''' pre-package in ColumnDataSources '''
galaxy_sb1 = ColumnDataSource(dict(x=galaxyfiles[0][0],y=galaxyfiles[0][1]))
galaxy_sb2 = ColumnDataSource(dict(x=galaxyfiles[1][0],y=galaxyfiles[1][1]))
galaxy_sb3 = ColumnDataSource(dict(x=galaxyfiles[2][0],y=galaxyfiles[2][1]))
galaxy_sb4 = ColumnDataSource(dict(x=galaxyfiles[3][0],y=galaxyfiles[3][1]))
galaxy_sb5 = ColumnDataSource(dict(x=galaxyfiles[4][0],y=galaxyfiles[4][1]))
galaxy_sb6 = ColumnDataSource(dict(x=galaxyfiles[5][0],y=galaxyfiles[5][1]))
galaxy_s0 = ColumnDataSource(dict(x=galaxyfiles[6][0],y=galaxyfiles[6][1]))
galaxy_sa = ColumnDataSource(dict(x=galaxyfiles[7][0],y=galaxyfiles[7][1]))
galaxy_sb = ColumnDataSource(dict(x=galaxyfiles[8][0],y=galaxyfiles[8][1]))
galaxy_sc = ColumnDataSource(dict(x=galaxyfiles[9][0],y=galaxyfiles[9][1]))
galaxy_bulge = ColumnDataSource(dict(x=galaxyfiles[10][0],y=galaxyfiles[10][1]))
galaxy_ellipticals = ColumnDataSource(dict(x=galaxyfiles[11][0],y=galaxyfiles[11][1]))
galaxy_lbg_all_flam_x = ColumnDataSource(dict(x=galaxyfiles[12][0],y=galaxyfiles[12][1]))

o5v = ColumnDataSource(dict(x=starfiles[0][0],y=starfiles[0][1]))
b0v = ColumnDataSource(dict(x=starfiles[1][0],y=starfiles[1][1]))
b57v = ColumnDataSource(dict(x=starfiles[2][0],y=starfiles[2][1]))
a0v = ColumnDataSource(dict(x=starfiles[3][0],y=starfiles[3][1]))
a5v = ColumnDataSource(dict(x=starfiles[4][0],y=starfiles[4][1]))
f0v = ColumnDataSource(dict(x=starfiles[5][0],y=starfiles[5][1]))
f5v = ColumnDataSource(dict(x=starfiles[6][0],y=starfiles[6][1]))
g0v = ColumnDataSource(dict(x=starfiles[7][0],y=starfiles[7][1]))
g5v = ColumnDataSource(dict(x=starfiles[8][0],y=starfiles[8][1]))
k0v = ColumnDataSource(dict(x=starfiles[9][0],y=starfiles[9][1]))
k5v = ColumnDataSource(dict(x=starfiles[10][0],y=starfiles[10][1]))
m0v = ColumnDataSource(dict(x=starfiles[11][0],y=starfiles[11][1]))
m5v = ColumnDataSource(dict(x=starfiles[12][0],y=starfiles[12][1]))
