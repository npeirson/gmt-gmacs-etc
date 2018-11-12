''' keys '''

'''
	example parameters dict object:
		my_params = dict(mode='snr',telescope_mode='first',object_type='a5v')

			)
'''

arguments = [['mode','type'],['num_mirrors''telescope_mode'],['wavelength','wavs'],['exposure_time','exp'],
			['object','object_type'],['filter_opt','filter'],['mag_sys','magnitude_system'],
			['mag','magnitude'],['redshift','z'],['seeing','astronomical_seeing'],['slit_width','slit_size'],
			['moon_days','new_moon'],['grating','grating_opt'],['read_noise','noise'],['bin_opt','binning'],
			['channel','channels'],['sss','verbose']]

mode = [['snr','signal_noise',None],['os','obs_spec','observed_spectrum'],
		['sky','sky_bg','sky_background'],['dichro','dichroic','dichroic_throughput'],
		['gt','grating','grating_throughput'],['ccd','ccd_qe','quantum_efficiency'],
		[r'atmo*_ext','atmospheric_extinction','extinction']]

telescope_mode = [['first','first_light',4],['full',None,7]]

object_type = [['SB1','SB2','SB3','SB4','SB5','SB6','S0','Sa','Sb','Sc','bulge','ellipticals','lbg_all_flam'],
				['o5v','b0v','b5v','a0v','a5v','f0v','g0v','g5v','k0v','k5v','m0v','m5v','f5v']]

filter_opt = ['u','g','r','i','z','_u','_b','_v','_r','_i']

mag_sys = ['vega','ab']

moon_days = [0,3,7,10,14]


keychain = [mode,telescope_mode,object_type,filter_opt]





stellar = np.concatenate(([filename[:-4] for filename in epaths.stellar_files],[i for i in range(len(epaths.stellar_files))]))
galactic = np.concatenate((epaths.galaxy_files,[(i+len(epaths.stellar_files)) for i in range(len(epaths.galaxy_files))]))