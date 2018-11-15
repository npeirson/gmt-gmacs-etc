#!/usr/bin/python3

import numpy as np
import paths as etpaths
import strings as stc


''' keys '''

arguments = [['mode','type'],['num_mirrors','telescope_mode'],['wavelength','wavs'],['exposure_time','exp'],
			['object_class','obj'],['filter_opt','filter'],['mag_sys','magnitude_system'],
			['mag','magnitude'],['redshift','z'],['seeing','astronomical_seeing'],['slit_width','slit_size'],
			['moon_days','new_moon'],['grating_opt','grating'],['readnoise','noise'],['bin_opt','binning'],
			['channel','channels'],['sss','silent']]
mode = [['snr','signal_noise',None],['os','obs_spec','observed_spectrum'],
		['sky','sky_bg','sky_background'],['dichro','dichroic','dichroic_throughput'],
		['gt','grating','grating_throughput'],['ccd','ccd_qe','quantum_efficiency'],
		[r'atmo*_ext','atmospheric_extinction','extinction']]
num_mirrors = [['first','first_light',4],['full',None,7]]
object_type = [['SB1','SB2','SB3','SB4','SB5','SB6','S0','Sa','Sb','Sc','bulge','ellipticals','lbg_all_flam'],
				['o5v','b0v','b5v','a0v','a5v','f0v','g0v','g5v','k0v','k5v','m0v','m5v','f5v']]
filter_opt = [['u','g','r','i','z','_u','_b','_v','_r','_i']]
mag_sys = [['vega'],['ab']]
moon_days = [[0,3,7,10,14]]
grating_opt = [['low','low_res','low_resolution',None],['high','high_res','high_resolution',True]]
bin_opt = [[1,2,3,4]]

keychain = [mode,num_mirrors,object_type,filter_opt,mag_sys,moon_days,grating_opt,moon_days,grating_opt,bin_opt]

stellar = [filename[:-4] for filename in etpaths.stellar_files]
galactic = etpaths.galaxy_files

valid_widgets = [name[7:] for name in stc.widget_names]


''' keychain '''
func_grating = ['grating','slit','init']
func_filter = ['filter','init'] # add plot_step if it's made dynamic
func_dichro = ['wavelength','init']
func_grating = ['wavelength','init']
func_ccd = ['wavelength','init']
func_mirror_loss = ['wavelength','init']
func_readnoise = ['slit','binning','seeing','wavelength','init']
func_atmo_ext = ['wavelength','init']

args_power = ['flux','area','time']
args_counts = np.concatenate((args_power,['wavelength']),axis=0)

args_percent = ['seeing','slit']

args_total_eff = ['wavelength']

args_signal = np.unique(np.concatenate((args_power,args_percent,args_total_eff),axis=0))

args_sky_flux = ['wavelength']
args_extension = ['seeing','slit']

args_counts_noise = np.unique(np.concatenate((args_sky_flux,args_extension,['time']),axis=0))
args_total_eff_noise = ['wavelength']

args_noise = np.unique(np.concatenate((args_counts_noise,args_total_eff_noise),axis=0))

args_snr = np.unique(np.concatenate((args_signal,args_noise),axis=0))

args_err = args_snr

args_os = np.unique(np.concatenate((args_signal,args_err),axis=0))

args_wavelength = ['wavelength','init']