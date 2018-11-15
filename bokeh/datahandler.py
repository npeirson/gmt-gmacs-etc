#!/usr/bin/python3

import os
import numpy as np
#import json
#import pandas as pd
import defaults as dfs
import paths as etpaths
#from bokeh.models import ColumnDataSource

# jsonify_data = lambda data,filename: json.dump(data,filename) # for an optimization in later versions

''' import data '''
# explicit... not very elegant, but pandas don't like loops
sb1_x,sb1_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[0]),usecols=(0,1),unpack=True)
sb2_x,sb2_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[1]),usecols=(0,1),unpack=True)
sb3_x,sb3_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[2]),usecols=(0,1),unpack=True)
sb4_x,sb4_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[3]),usecols=(0,1),unpack=True)
sb5_x,sb5_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[4]),usecols=(0,1),unpack=True)
sb6_x,sb6_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[5]),usecols=(0,1),unpack=True)
s0_x,s0_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[6]),usecols=(0,1),unpack=True)
sa_x,sa_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[7]),usecols=(0,1),unpack=True)
sb_x,sb_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[8]),usecols=(0,1),unpack=True)
sc_x,sc_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[9]),usecols=(0,1),unpack=True)
bulge_x,bulge_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[10]),usecols=(0,1),unpack=True)
ellipticals_x,ellipticals_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[11]),usecols=(0,1),unpack=True)
lbg_all_flam_x,lbg_all_flam_y = np.loadtxt(os.path.join(etpaths.galaxy_path,etpaths.galaxy_files[12]),usecols=(0,1),unpack=True)

o5v_x,o5v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[0]),usecols=(0,1),unpack=True)
b0v_x,b0v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[1]),usecols=(0,1),unpack=True)
b57v_x,b57v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[2]),usecols=(0,1),unpack=True)
a0v_x,a0v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[3]),usecols=(0,1),unpack=True)
a5v_x,a5v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[4]),usecols=(0,1),unpack=True)
f0v_x,f0v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[5]),usecols=(0,1),unpack=True)
f5v_x,f5v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[12]),usecols=(0,1),unpack=True)
g0v_x,g0v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[6]),usecols=(0,1),unpack=True)
g5v_x,g5v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[7]),usecols=(0,1),unpack=True)
k0v_x,k0v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[8]),usecols=(0,1),unpack=True)
k5v_x,k5v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[9]),usecols=(0,1),unpack=True)
m0v_x,m0v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[10]),usecols=(0,1),unpack=True)
m5v_x,m5v_y = np.loadtxt(os.path.join(etpaths.stellar_path,etpaths.stellar_files[11]),usecols=(0,1),unpack=True)

filter_besb_x,filter_besb_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[0]),usecols=(0,1),unpack=True,delimiter=',')
filter_besi_x,filter_besi_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[1]),usecols=(0,1),unpack=True,delimiter=',')
filter_besr_x,filter_besr_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[2]),usecols=(0,1),unpack=True,delimiter=',')
filter_besu_x,filter_besu_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[3]),usecols=(0,1),unpack=True,delimiter=',')
filter_besv_x,filter_besv_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[4]),usecols=(0,1),unpack=True,delimiter=',')
filter_u_x,filter_u_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[5]),usecols=(0,1),unpack=True,delimiter=',')
filter_g_x,filter_g_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[6]),usecols=(0,1),unpack=True,delimiter=',')
filter_r_x,filter_r_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[7]),usecols=(0,1),unpack=True,delimiter=',')
filter_i_x,filter_i_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[8]),usecols=(0,1),unpack=True,delimiter=',')
filter_z_x,filter_z_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[9]),usecols=(0,1),unpack=True,delimiter=',')
filter_photonux_x,filter_photonux_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[10]),usecols=(0,1),unpack=True)
filter_photonb_x,filter_photonb_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[11]),usecols=(0,1),unpack=True)
filter_photonv_x,filter_photonv_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[12]),usecols=(0,1),unpack=True)
filter_photonr_x,filter_photonr_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[13]),usecols=(0,1),unpack=True)
filter_photoni_x,filter_photoni_y = np.loadtxt(os.path.join(etpaths.filter_path,etpaths.filter_files[14]),usecols=(0,1),unpack=True)


skyfile_00d_x,skyfile_00d_y = np.loadtxt(os.path.join(etpaths.skyfiles_path,etpaths.skyfiles[0]),usecols=(0,1),unpack=True)
skyfile_03d_x,skyfile_03d_y = np.loadtxt(os.path.join(etpaths.skyfiles_path,etpaths.skyfiles[1]),usecols=(0,1),unpack=True)
skyfile_07d_x,skyfile_07d_y = np.loadtxt(os.path.join(etpaths.skyfiles_path,etpaths.skyfiles[2]),usecols=(0,1),unpack=True)
skyfile_10d_x,skyfile_10d_y = np.loadtxt(os.path.join(etpaths.skyfiles_path,etpaths.skyfiles[3]),usecols=(0,1),unpack=True)
skyfile_14d_x,skyfile_14d_y = np.loadtxt(os.path.join(etpaths.skyfiles_path,etpaths.skyfiles[4]),usecols=(0,1),unpack=True)

grating_red_x,grating_red_y = np.loadtxt(etpaths.grating_path[0],usecols=(0,1),unpack=True)
grating_blue_x,grating_blue_y = np.loadtxt(etpaths.grating_path[1],usecols=(0,1),unpack=True)
ccd_red_x,ccd_red_y = np.loadtxt(etpaths.ccd_path[0],usecols=(0,1),unpack=True)
ccd_blue_x,ccd_blue_y = np.loadtxt(etpaths.ccd_path[1],usecols=(0,3),unpack=True)

_dichroic_x,dichroic_y1,dichroic_y2 = np.loadtxt(etpaths.dichroic_path,usecols=(0,1,2),unpack=True)
dichroic_x = _dichroic_x * 10
atmo_ext_x,atmo_ext_y = np.loadtxt(etpaths.atmo_ext_path,usecols=(0,1),unpack=True,delimiter=',')

mirror_file = np.fliplr(np.loadtxt(etpaths.mirror_file,usecols=(0,1),unpack=True,delimiter=' '))
vega_file = np.loadtxt(etpaths.vega_file,usecols=(0,1),unpack=True,delimiter=',')

''' coalesced data '''

galaxyfiles = [[sb1_x,sb1_y],[sb2_x,sb2_y],[sb3_x,sb3_y],[sb4_x,sb4_y],[sb5_x,sb5_y],
				[sb6_x,sb6_y],[s0_x,s0_y],[sa_x,sa_y],[sb_x,sb_y],[sc_x,sc_y],
				[bulge_x,bulge_y],[ellipticals_x,ellipticals_y],[lbg_all_flam_x,lbg_all_flam_y]]
starfiles = [[o5v_x,o5v_y],[b0v_x,b0v_y],[b57v_x,b57v_y],
			[a0v_x,a0v_y],[a5v_x,a5v_y],[f0v_x,f0v_y],[f5v_x,f5v_y],
			[g0v_x,g0v_y],[g5v_x,g5v_y],[k0v_x,k0v_y],[k5v_x,k5v_y],
			[m0v_x,m0v_y],[m5v_x,m5v_y]]
filterfiles = [[filter_besb_x,filter_besb_y],[filter_besi_x,filter_besi_y],
				[filter_besr_x,filter_besr_y],[filter_besu_x,filter_besu_y],
				[filter_besv_x,filter_besv_y],[filter_photonux_x,filter_photonux_y],
				[filter_photonb_x,filter_photonb_y],[filter_photonv_x,filter_photonv_y],
				[filter_photonr_x,filter_photonr_y],[filter_photoni_x,filter_photoni_y],
				[filter_u_x,filter_u_y],[filter_g_x,filter_g_y],[filter_r_x,filter_r_y],
				[filter_i_x,filter_i_y],[filter_z_x,filter_z_y]]
skyfiles = [[skyfile_00d_x,skyfile_00d_y],[skyfile_03d_x,skyfile_03d_y],[skyfile_07d_x,skyfile_07d_y],
			[skyfile_10d_x,skyfile_10d_y],[skyfile_14d_x,skyfile_14d_y]]
dichroic_x,dichroic_y1,dichroic_y2 = dichroic_x,dichroic_y1,dichroic_y2
grating1,grating2 = [grating_blue_x,grating_blue_y],[grating_red_x,grating_red_y]
ccd1,ccd2 = [ccd_blue_x,ccd_blue_y],[ccd_red_x,ccd_red_y]
atmo_ext_x,atmo_ext_y = atmo_ext_x,atmo_ext_y
mirror_file_x,mirror_file_y = mirror_file[0]*10,mirror_file[1]