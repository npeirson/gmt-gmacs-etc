# this uses a caller singleton and enumerated caller-association groups to update only necessary data.
# the else statement catches errors and wavelength changes.


import numpy as np
import values as edl
# coating efficiencies moved into values file, edl.coating_eff_red

def type_router(_obj):
	try:
		if isinstance(_obj,str):
			return int(_obj.lower())
		elif isinstance(_obj,int):
			return _obj
		else:
			raise ValueError("{} Invalid value: {}".format(string_prefix,_obj))
	except:
		raise TypeError("{} Invalid type: {}".format(string_prefix,type(_obj)))

class simulate:
	string_prefix = '[ etc ] :'
	fmirror = interpolate.interp1d(mirror_file_x,mirror_file_y, kind='cubic')
	mirror = fmirror(wavelength)

	def __init__(self,telescope_mode='first',wavelength=np.arange(3200,10360,edl.dld[0]/3.),exposure_time=3600,object_type='a5v',
		filter_index=3,mag_sys_opt='ab',magnitude=25,redshift=0,seeing=0.5,slit_size=0.5,moon_days=0,grating_opt=0,noise=False,
		bin_option=edl.bin_options_int[edl.bin_options_default_index],channel='both'):
		
		object_type = type_router(object_type)
		grating_opt = type_router(grating_opt)
		telescope_mode = type_router(telescope_mode)
		self.plot_step = wavelength[2] - wavelength[1]


	''' tier 0 '''

	def refresh(self,caller,tabs=tabs): # default cb_obj as immutable
		# direct the changed values
		if caller in edl.signal_components:
			recalculate_signal(caller)
		if caller in edl.noise_components:
			recalculate_noise(caller)
		else:
			recalculate_signal(caller)
			recalculate_noise(caller)
		# recalculate only necessary stuff
		if (tabs.active == 0): # snr
			recalculate_snr(caller)
		elif (tabs.active == 1): # obs spec no noise
			pass # just send signal back :)
		elif (tabs.active == 2): # obs spec noise
			recalculate_error(caller)
		elif (tabs.active == 3): #obs sky background
			pass # send noise
		return the_stuff


	def recalculate_snr(self,caller):
		if (self.channel == 'blue') or (self.channel == 'both'):
			self.snr_blue = np.divide(self.blue_signal,np.sqrt(self.blue_signal + self.blue_noise + np.square(self.readnoise)))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.snr_red = np.divide(self.red_signal,np.sqrt(self.red_signal + self.red_noise + np.square(self.readnoise)))


	def recalculate_error(self,caller):
		if (self.self.channel == 'blue') or (self.channel == 'both'):
			sigma_blue = np.sqrt(self.blue_signal + self.blue_noise + np.square(self.readnoise))
			self.error_blue = np.random.normal(loc=0, scale=sigma_blue,size=len(self.snr_blue))
		if (self.channel == 'red') or (self.channel == 'both'):
			sigma_red = np.sqrt(self.red_signal + self.red_noise + np.square(self.readnoise))
			self.error_red = np.random.normal(loc=0, scale=sigma_red,size=len(self.snr_red))


	''' tier 1 '''

	def recalculate_signal(self,caller):
		if caller in edl.counts_components:
			pass # recalc counts stuff
		elif caller in edl.percent_components:
			pass # recalc percent stuff
		elif caller in edl.efficiency_components:
			pass # recalc efficiency stuff
		else:
			pass # recalc all stuff

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.blue_signal = np.multiply((self.counts * self.percent), self.blue_total_eff)
		if (self.channel == 'red') or (self.channel == 'both'):
			self.red_signal = np.multiply((self.counts * self.percent), self.red_total_eff)

	def recalculate_noise(self,caller):
		if caller in edl.counts_noise_components:
			pass # recalc counts noise stuff
		elif caller in edl.efficiency_noise_components:
			pass # recalc total efficiency noise stuff

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.blue_noise = np.multiply(self.counts_noise,self.blue_total_eff_noise)
		if (channel == 'red') or (self.channel == 'both'):
			self.red_noise = np.multiply(self.counts_noise,self.red_total_eff_noise)

	''' tier 2 '''

	def recalculate_counts(self,caller):
		flux = self.flux_y
		self.power = flux * self.area * self.exposure_time * self.plot_step
		self.counts = np.divide(np.divide(self.power,np.divide((const.h.value * const.c.value),self.wavelength)),1e10)


	def recalculate_counts_noise(self,caller):
		recalculate_sky_flux(caller)
		recalculate_extension(caller)
		self.counts_noise = np.multiply(np.multiply(self.sky_flux,self.extension),(self.area*self.exposure_time*self.plot_step))


	def recalculate_percent(self,caller):
		recalculate_seeing(caller)
		# this is a forwarding function to help me keep my head-thoughts straight


	def recalculate_efficiency(self,caller):
		if (self.channel == 'blue') or (self.channel == 'both'):
			self.blue_total_eff = np.multiply(np.multiply(self.blue_dichro,self.blue_grating),np.multiply((blue_ccd * (edl.coating_eff_blue * self.extinction)),np.square(mirror)))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.red_total_eff = np.multiply(np.multiply(self.red_dichro,self.red_grating),np.multiply((self.red_ccd * (edl.coating_eff_red * self.extinction)),np.square(mirror)))


	def recalculalte_efficiency_noise(self,caller):
		if (self.channel == 'blue') or (self.channel == 'both'):
			self.blue_total_eff_noise = np.multiply(np.multiply(self.blue_dichro,self.blue_grating),(self.blue_ccd * np.square(self.mirror) * edl.coating_eff_blue))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.red_total_eff_noise = np.multiply(np.multiply(self.red_dichro,self.red_grating),(self.red_ccd * np.square(self.mirror) * edl.coating_eff_red))


	def recalculate_flux(self,caller):
		# heal identicalities
		lambda_A[0] = lambda_A[0] + plot_step
		lambda_A[-1] = lambda_A[-1] - plot_step

		ftrans = interpolate.interp1d(selected_filter[0],selected_filter[1], kind='cubic')
		trans = ftrans(lambda_A) #spectres(lambda_A,selected_filter[0],selected_filter[1])

		_extinction = spectres(lambda_A,atmo_ext_x,atmo_ext_y)

		flux = flux_A * 1e10
		_lambda = lambda_A / 1e10

		num_zeros = 0
		for lux in flux:
			if (lux is None) or (lux is 0):
				num_zeros += 1
				lux = 0

		if (num_zeros >= (flux.shape[0]/5)):
			if (num_zeros == flux.shape[0]):
				print('No flux in this bandpass!')
				output_flux = np.zeros(wavelength.shape[0])
			else:
				percent_zeros = (num_zeros / flux.shape[0]) * 100
				print('{}% of this bandpass has zero flux'.format(percent_zeros))

		if (mag_sys_opt == 'vega'):
			flux_vega = spectres(wavelength,vega[0],vega[1]) * 1e10 # probably not right Need to read in vega file, not defined right now
			print(vega[0])
			mag_model = -2.5 * np.log10(np.divide(math.fsum(flux * _extinction * _lambda * trans),math.fsum(flux_vega * trans * _lambda * _extinction))) + 0.03
		elif (mag_sys_opt == 'ab'):
			mag_model = -48.6 - 2.5 * np.log10(math.fsum(flux * trans * _extinction *_lambda) / math.fsum(trans * _lambda * _extinction * (const.c.value/np.square(_lambda))))
		else:
			print('Invalid mag_sys_opt!')

		del_mag = mag - mag_model
		output_lambda = object_x
		output_flux = np.multiply(object_y,10 ** np.negative(del_mag/2.5))
		old_res = output_lambda[2] - output_lambda[1]
		if (old_res < plot_step):
			self.flux_y = spectres(self.wavelength,object_x,(output_flux*1e-03)) # ergs s-1 cm-2 A-1 to J s-1 m-2 A-1
		else:
			self.flux_y = spectres(self.wavelength,object_x,(output_flux*1e-03))


	''' tier 3 '''

	def change_grating_opt(self,caller):
		if self.grating_opt in edl.grating_opt_keys[:2]: # low resolution
			self.delta_lambda = edl.dld[0] * self.slit_size / 0.7
		elif self.grating_opt in edl.grating_opt_keys[3:]:
			self.delta_lambda = edl.dld[1] * self.slit_size / 0.7
		else:
			raise ValueError("{} Invalid grating_opt: {}".format(string_prefix,self.grating_opt))


	def change_object_type(self,caller):
		if self.object_type in edl.stellar_keys:
			index_of = [i for i,name in enumerate(stellar_keys) if object_type in name][0]
			self.object_type = edl.starfiles[index_of]
		elif self.object_type in edl.galactic_keys:
			index_of = [i for i,name in enumerate(galactic_keys) if object_type in name][0]
			self.object_type = edl.galaxyfiles[index_of]
		else:
			raise ValueError("{} Invalid object type: {}".format(string_prefix,self.object_type))

		self.object_x = object_type[0] * (1+redshift)
		self.object_y = object_type[1]
		self.flux_A = spectres(lambda_A,object_x,object_y)


	def change_moon_days(self,caller):
		if moon_days in moon_days_keys:
			self.sky_background = skyfiles[(int(np.where(np.asarray(moon_days_keys)==moon_days)[0]))]
		else:
			raise ValueError('{} Invalid number of days since new moon: {}'.format(string_prefix,self.moon_days))
		recalculate_counts_noise(caller=caller)


	def change_telescope_mode(self,caller):
		if self.telescope_mode in edl.telescope_mode_keys[:3]:
			self.area = edl.area[0]
		elif self.telescope_mode in edl.telescope_mode_keys[4:]:
			self.area = edl.area[1]
		else:
			raise ValueError('{} Invalid telescope mode: {}'.format(string_prefix,self.telescope_mode))


	def change_filter(self,caller):
		selected_filter = dh.filterfiles[self.filter_index]
		filter_min = min(selected_filter[0])
		filter_max = max(selected_filter[0])

		if (filter_min > self.wavelength[0]):
			lambda_min = filter_min
		elif (filter_min == self.wavelength[0]):
			filter_min = selected_filter[int(np.where(selected_filter[0] > self.wavelength[0])[0])]
		else:
			lambda_min = self.wavelength[0]

		if (filter_max < self.wavelength[-1]):
			lambda_max = filter_max
		elif (filter_max == self.wavelength[-1]):
			filter_max = selected_filter[int(np.where(selected_filter[0] < self.wavelength[-1])[-1])]
		else:
			lambda_max = self.wavelength[-1]

		self.lambda_A = np.arange(lambda_min,lambda_max,self.plot_step)


	def recalculate_seeing(self,caller):
		_sigma = self.seeing / gaussian_sigma_to_fwhm
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		self.percent_u,self.percent_err_u = integrate.quad(funx,(-slit_size/2),(slit_size/2))
		self.percent_l,self.percent_err_l = integrate.quad(funx,(-seeing/2),(seeing/2))
		self.percent = self.percent_u * self.percent_l # can use error if you add it later...
		self.extension = self.seeing * self.slit_size


	def recalculate_sky_flux(self,caller):
		sky_x = self.sky_background[0] * 1e4
		sky_y = self.sky_background[1] / 1e4
		old_res = sky_x[2] - sky_x[1]
		_sigma = self.delta_lambda / gaussian_sigma_to_fwhm
		_x = np.arange((-5*_sigma),(5*_sigma),old_res)
		degrade = funx(_x)/np.trapz(funx(_x))
		sky_y = convolve_fft(sky_y,degrade)
		self.sky_flux = spectres(self.wavelength,sky_x,sky_y)


	def recalculate_dichroic(self,caller):
		if (self.channel is 'blue') or (self.channel is 'both'):
			fblue_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y1, kind='cubic')
			self.blue_dichro = fblue_dichro(self.wavelength)
		if (self.channel is 'red') or (self.channel is 'both'):
			fred_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y2, kind='cubic')
			self.red_dichro = fred_dichro(self.wavelength)


	def recalculate_grating(self,caller):
		if (self.channel is 'blue') or (self.channel is 'both'):
			self.blue_grating = spectres(self.wavelength,(dh.grating1[0]*10),dh.grating1[1])
		if (self.channel is 'red') or (self.channel is 'both'):
			self.red_grating = spectres(self.wavelength,(dh.grating2[0]*10),dh.grating2[1])


	def recalculate_ccd(self,caller):
		if (channel is 'blue') or (channel is 'both'):
			fblue_ccd = interpolate.interp1d((dh.ccd1[0]*10),dh.ccd1[1], kind='cubic')
			self.blue_ccd = fblue_ccd(self.wavelength)
		if (channel is 'red') or (channel is 'both'):
			fred_ccd = interpolate.interp1d((dh.ccd2[0]*10),dh.ccd2[1], kind='cubic')
			self.red_ccd = fred_ccd(wavelength)

	def recalculate_readnoise(self,caller):
		spectral_resolution = math.ceil((self.slit_size/(0.7/12))/2)*2 #px (ceil()/2)*2 to round up to next even integer
		spatial_resolution = math.ceil((self.seeing/(0.7/12))/2)*2 #px (ceil()/2)*2 to round up to next even integer 
		extent = seeing * slit_size
		npix = extent/(0.7/12)**2
		
		rn = edl.rn_default
		if (self.bin_size > 0) and (self.bin_size < 5):
			print('[ info ] : Pixel binning: ({}x{})'.format(self.bin_size,self.bin_size))
			readnoise = math.ceil(rn * spectral_resolution * spatial_resolution / (self.bin_size**2))
			print('[ info ] : Extent: {} arcsec^2\n[ info ] : num pixels/resel: {} px\n[ info ] : spectral resolution: {} px\n[ info ] : spatial resolution: {} px'.format(extent,int(math.ceil(npix)),spectral_resolution,spatial_resolution))
		else:
			raise ValueError('{} Invalid pixel binning option ({})'.format(string_prefix,self.bin_size))

	def recalculate_atmospheric_extinction(self,caller):
		self.extinction = spectres(self.wavelength,dh.atmo_ext_x,dh.atmo_ext_y)