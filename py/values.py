'''
	values file
		for easy changing and stuffs

		GMT info here: https://www.gmto.org/overview/quick-facts/
		GMACS info here: http://instrumentation.tamu.edu/gmacs.html

'''
import os

''' constants '''
coating_efficiency = 0.8*0.98**14
area = [222,368] #
dld = [3.73,1.4]
rn_default = 2
bin_options_default_index = 3 # 1th index is 2nd cell, ergo default is 2x2 binning
coating_eff_red = 0.62
coating_eff_blue = 0.60
bin_options_int = [1,2,3,4]
bin_options_str = ['1x1','2x2','3x3','4x4']

''' paths '''
# get path to values.py
absFilePath = os.path.abspath(__file__)
# get directory
fileDir = os.path.dirname(os.path.abspath(__file__))
# get parent directory
parentDir = os.path.dirname(fileDir)

galaxy_path = parentDir + '/core/kinney/'
stellar_path = parentDir + '/core/pickle/'
skyfiles_path = parentDir + '/core/skybackground/'
filter_path = parentDir + '/core/filters/'
atmo_ext_path = parentDir + '/core/atmo_extinction.dat'
dichroic_path = parentDir + '/core/dichroic5577.txt'
grating_path = [parentDir + '/core/grating_red_low_315-1100.txt',parentDir + '/core/grating_blue_low_315-1100.txt']
ccd_path = [parentDir + '/core/e2v-astro-multi-2-DD.txt',parentDir + '/core/e2v_blue.txt']
vega_file = parentDir + '/core/alpha_lyr_stis_005.ascii'
mirror_file = parentDir + '/core/mirror_300-1200nm.txt'
galaxy_files = ['SB1','SB2','SB3','SB4','SB5','SB6','S0','Sa','Sb','Sc','bulge','ellipticals','lbg_all_flam']
stellar_files = ['o5v.dat','b0v.dat','b57v.dat','a0v.dat','a5v.dat','f0v.dat','g0v.dat','g5v.dat','k0v.dat','k5v.dat','m0v.dat','m5v.dat']
skyfiles = ['00d_315-1200nm.csv','03d_315-1200nm.csv','07d_315-1200nm.csv','10d_315-1200nm.csv','14d_315-1200nm.csv']
#filter_files = ['BesB.dat','BesI.dat','BesR.dat','BesU.dat','BesV.dat','g.dat','i.dat','photonb.dat','photonI.dat','photonR.dat','photonUX.dat','photonV.dat','r.dat','u.dat','z.dat']
filter_files = ['photonUX.dat','photonB.dat','photonV.dat','photonR.dat','photonI.dat','u.dat','g.dat','r.dat','i.dat','z.dat']

galaxyfiles = [[sb1_x,sb1_y],[sb2_x,sb2_y],[sb3_x,sb3_y],[sb4_x,sb4_y],[sb5_x,sb5_y],
				[sb6_x,sb6_y],[s0_x,s0_y],[sa_x,sa_y],[sb_x,sb_y],[sc_x,sc_y],
				[bulge_x,bulge_y],[ellipticals_x,ellipticals_y],[lbg_all_flam_x,lbg_all_flam_y]]
starfiles = [[star_o5v_x,star_o5v_y],[star_b0v_x,star_b0v_y],[star_b57v_x,star_b57v_y],
			[star_a0v_x,star_a0v_y],[star_a5v_x,star_a5v_y],[star_f0v_x,star_f0v_y],
			[star_g0v_x,star_g0v_y],[star_g5v_x,star_g5v_y],[star_k0v_x,star_k5v_y],
			[star_m0v_x,star_m0v_y],[star_m5v_x,star_m5v_y]]
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

''' keychain '''
# argument keys (probably not right, yet)
keys = ['obj','atmo_ext','wavelength',
		'filter_opt','magnitude','mag_sys_opt',
		'grating_opt','redshift','exposure_time',
		'seeing','slit_size','moon_days','plot_channel',
		'telescope_mode','binx','sss']

# passage keys, wavelength assumed
stellar_keys = [filename[:-4] for filename in stellar_files]
galactic_keys = edl.galaxy_files
filter_keys = [filename[:-4] for filename in filter_files]
grating_opt_keys = ['low',0,1.4,'high',1,3.73]
moon_days_keys = [0,3,7,10,14]
telescope_mode_keys = [0,4,'first','first light',1,7,'full','full size']

power_keys = ['flux','area','exposure_time']
percent_keys = ['seeing']
total_eff_keys = ['grating','dichro','ccd','atmo_ext']
total_eff_noise_keys = []
obs_sky_background_keys = ['noise']
dichroic_keys = []
grating_keys = ['grating_opt','slit_size']
ccd_keys = []
atmo_ext_keys = []
extension_keys = ['seeing','slit_size']
mirror_keys = []

# `np.unique` makes it safe to list key-dependencies in their entirety, for mental health
counts_noise_keys = np.unique(np.concatenate(['moon_days','telescope_mode','exposure_time'],extension_keys))
flux_keys = np.unique(np.concatenate(['grating_opt','filter'],filter_keys))
counts_keys = np.unique(np.concatenate(['telescope_mode'],power_keys))
signal_keys = np.unique(np.concatenate(total_eff_keys,counts_keys,percent_keys))
noise_keys = np.unique(np.concatenate(['moon_days'],counts_noise_keys,total_eff_noise_keys))
error_keys = np.unique(np.concatenate(signal_keys,noise_keys))
readnoise_keys = np.unique(np.concatenate(['binning'],extension_keys))

snr_keys = np.unique(np.concatenate(signal_keys,noise_keys,readnoise_keys))
obs_spec_noise_keys = np.unique(np.concatenate(signal_keys,error_keys))
obs_spec_nonoise_keys = np.unique(np.concatenate(signal_keys)))