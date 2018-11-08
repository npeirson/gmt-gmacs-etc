'''
	values file
		for easy changing and stuffs

		GMT info here: https://www.gmto.org/overview/quick-facts/
		GMACS info here: http://instrumentation.tamu.edu/gmacs.html

'''
import os
import numpy as np

''' constants '''
string_prefix = '[ etc ] :'
coating_efficiency = 0.8*0.98**14
area = [222,368]
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
stellar_files = ['o5v.dat','b0v.dat','b57v.dat','a0v.dat','a5v.dat','f0v.dat','g0v.dat','g5v.dat','k0v.dat','k5v.dat','m0v.dat','m5v.dat','f5v.dat']
skyfiles = ['00d_315-1200nm.csv','03d_315-1200nm.csv','07d_315-1200nm.csv','10d_315-1200nm.csv','14d_315-1200nm.csv']
#filter_files = ['BesB.dat','BesI.dat','BesR.dat','BesU.dat','BesV.dat','g.dat','i.dat','photonb.dat','photonI.dat','photonR.dat','photonUX.dat','photonV.dat','r.dat','u.dat','z.dat']
filter_files = ['photonUX.dat','photonB.dat','photonV.dat','photonR.dat','photonI.dat','u.dat','g.dat','r.dat','i.dat','z.dat']


''' keychain '''
# argument keys (probably not right, yet)
keys = ['obj','atmo_ext','wavelength',
		'filter_opt','magnitude','mag_sys_opt',
		'grating_opt','redshift','exposure_time',
		'seeing','slit_size','moon_days','plot_channel',
		'telescope_mode','binx','sss']

# passage keys, wavelength assumed
stellar_keys = [filename[:-4] for filename in stellar_files]
galactic_keys = galaxy_files
filter_keys = [filename[:-4] for filename in filter_files]
object_type_keys = np.concatenate((stellar_keys,galactic_keys))
grating_opt_keys = ['low',0,1.4,'high',1,3.73]
moon_days_keys = [0,3,7,10,14]
telescope_mode_keys = [0,4,'first','first light',1,7,'full','full size']

mag_sys_opt_keys = ['flux'] # encompasses much more than just flux, though
power_keys = ['flux','area','exposure_time']
percent_keys = ['seeing']
total_eff_keys = ['grating','dichro','ccd','atmo_ext']
obs_sky_background_keys = ['noise']
grating_keys = ['grating_opt','slit_size']
extension_keys = ['seeing','slit_size']

# `np.unique` makes it safe to list key-dependencies in their entirety, for mental health
counts_noise_keys = np.unique(np.concatenate((['moon_days','telescope_mode','exposure_time'],extension_keys)))
flux_keys = np.unique(np.concatenate((['grating_opt','filter'],filter_keys)))
counts_keys = np.unique(np.concatenate((['telescope_mode'],power_keys)))
signal_keys = np.unique(np.concatenate((total_eff_keys,counts_keys,percent_keys)))
noise_keys = np.unique(np.concatenate((['moon_days'],counts_noise_keys)))
error_keys = np.unique(np.concatenate((signal_keys,noise_keys)))
readnoise_keys = np.unique(np.concatenate((['binning'],extension_keys)))

snr_keys = np.unique(np.concatenate((signal_keys,noise_keys,readnoise_keys)))
obs_spec_noise_keys = np.unique(np.concatenate((signal_keys,error_keys)))
obs_spec_nonoise_keys = signal_keys


'''
initial_values = dict(gly_snr_red=gly_snr_red,gly_os_noise_red=gly_os_noise_red,gly_os_nonoise_red=gly_os_nonoise_red,
	gly_sky_red=gly_sky_red,gly_dichroic_red=gly_dichroic_red,gly_grating_red=gly_grating_red,gly_ccd_red=gly_ccd_red,
	gly_snr_blue=gly_snr_blue,gly_os_noise_blue=gly_os_noise_blue,gly_os_nonoise_blue=gly_os_nonoise_blue,
	gly_sky_blue=gly_sky_blue,gly_dichroic_blue=gly_dichroic_blue,gly_grating_blue=gly_grating_blue,gly_ccd_blue=gly_ccd_blue,
	tabs=tabs,telescope_mode='first',wavelength=dfs.default_wavelength,exposure_time=3600,object_type='a5v',
	filter_index=3,mag_sys_opt='ab',magnitude=25,redshift=0,seeing=0.5,slit_size=0.5,moon_days=0,grating_opt=0,
	noise=False,bin_option=edl.bin_options_int[edl.bin_options_default_index],channel='both',sss=True)
]
'''