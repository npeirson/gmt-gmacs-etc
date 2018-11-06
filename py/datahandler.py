import os
import json
import pandas as pd
import numpy as np
import config as cfg
import defaults as dfs
import values as edl


jsonify_data = lambda data,filename: json.dump(data,filename)

''' data handling '''
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

star_o5v_x,star_o5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[0]),usecols=(0,1),unpack=True)
star_b0v_x,star_b0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[1]),usecols=(0,1),unpack=True)
star_b57v_x,star_b57v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[2]),usecols=(0,1),unpack=True)
star_a0v_x,star_a0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[3]),usecols=(0,1),unpack=True)
star_a5v_x,star_a5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[4]),usecols=(0,1),unpack=True)
star_f0v_x,star_f0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[5]),usecols=(0,1),unpack=True)
star_g0v_x,star_g0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[6]),usecols=(0,1),unpack=True)
star_g5v_x,star_g5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[7]),usecols=(0,1),unpack=True)
star_k0v_x,star_k0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[8]),usecols=(0,1),unpack=True)
star_k5v_x,star_k5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[9]),usecols=(0,1),unpack=True)
star_m0v_x,star_m0v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[10]),usecols=(0,1),unpack=True)
star_m5v_x,star_m5v_y = np.loadtxt(os.path.join(edl.stellar_path,edl.stellar_files[11]),usecols=(0,1),unpack=True)

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