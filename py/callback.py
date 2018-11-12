# handles various input forms
class simulate:
	"""
	Class for simulating observatations using GMACS, designed to serve as the back-end to a static bokeh site.

	Parameters
	----------
		todo


	"""
	def __init__(self,bokeh_inserts=bokeh_inserts,initial_values=initial_values,**kwargs):
		# required args, maybe add range limit too
		'''
		self.object_type = self.type_router(initial_values['object_type'],edl.object_type_keys)
		if (self.object_type >= len(edl.object_type_keys)):
			self.object_type = self.object_type - len(edl.object_type_keys)
		'''
		self.old_mso = ''
		self.initialize_singletons = True

		try:
			if isinstance(initial_values['sss'],bool):
				self.sss = initial_values['sss']
			else:
				raise TypeError('{} Subsiste sermonem statim must be type bool, not type {}'.format(self.string_prefix,type(initial_values['sss'])))
		except:
			self.sss = True

		try:
			if isinstance(initial_values['channel'],str):
				if initial_values['channel'] in edl.channel_keys:
					self.channel = initial_values['channel']
				else:
					self.channel = 'both'
			else:
				raise ValueError("{} Invalid initial channel option: {}".format(edl.string_prefix,initial_values['channel']))
		except:
			print('uh oh...')
			self.channel = 'both' # improve this

		self.bin_size = 2 # this too


		try:
			if isinstance(initial_values['exposure_time'],int) or isinstance(initial_values['exposure_time'],float):
				self.exposure_time = initial_values['exposure_time']
			else:
				raise ValueError('{} Invalid exposure time: {}'.format(edl.string_prefix,self.exposure_time))
		except: # fix these error messages, they're not really... right... technically.
			raise TypeError('{} Invalid exposure time {}'.format(edl.string_prefix,initial_values['exposure_time']))

		self.object_type = initial_values['object_type']

		self.grating_opt = self.type_router(initial_values['grating_opt'],edl.grating_opt_keys)
		if (self.grating_opt >= 3):
			self.grating_opt = 1
		else:
			self.grating_opt = 0

		self.telescope_mode = self.type_router(initial_values['telescope_mode'],edl.telescope_mode_keys)
		if (self.telescope_mode >= 3):
			edl.telescope_mode = 'full'
		else:
			edl.telescope_mode = 'first'

		try:
			if isinstance(initial_values['mag_sys_opt'],str):
				if (initial_values['mag_sys_opt'].lower() == 'vega') or (initial_values['mag_sys_opt'].lower() == 'ab'):
					self.mag_sys_opt = initial_values['mag_sys_opt']
				else:
					raise ValueError('{} Invalid magnitude system option (mag_sys_opt) {}'.format(edl.string_prefix,initial_values['mag_sys_opt'].lower()))
			elif isinstance(initial_values['mag_sys_opt'],int) or isinstance(initial_values['mag_sys_opt'],float):
				if (int(initial_values['mag_sys_opt']) == 0):
					self.mag_sys_opt = 'vega'
				elif (int(initial_values['mag_sys_opt']) == 1):
					self.mag_sys_opt = 'ab'
				else:
					raise ValueError('{} Invalid magnitude system option (mag_sys_opt) {}'.format(edl.string_prefix,str(initial_values['mag_sys_opt'])))
			if not self.sss:
				print('{} Magnitude system set to {}'.format(edl.string_prefix,self.mag_sys_opt))
		except:
			self.mag_sys_opt = 'ab' # default
			print('{} Defaulting to magnitude system option {}'.format(edl.string_prefix,self.mag_sys_opt))

		self.filter_index = self.type_router(initial_values['filter_index'],edl.filter_keys)
		self.magnitude = self.num_router(initial_values['magnitude'])
		self.seeing = self.num_router(initial_values['seeing'])
		self.slit_size = self.num_router(initial_values['slit_size'])
		self.redshift = self.num_router(initial_values['redshift'])
		self.moon_days = self.num_router(initial_values['moon_days'])

		if isinstance(initial_values['wavelength'],np.ndarray):
			self.wavelength = initial_values['wavelength']
			if (self.wavelength[0] >= dfs.default_limits_wavelength[0]) and (self.wavelength[-1] <= dfs.default_limits_wavelength[1]):
				self.plot_step = self.wavelength[2] - self.wavelength[1]
			else:
				raise ValueError('{} Invalid wavelength extrema ({}--{}). Must be within {}--{}'.format(edl.string_prefix,wavelength[0],wavelength[-1],dfs.default_limits_wavelength[0],dfs.default_limits_wavelength[1]))
		else:
			raise TypeError("{} Invalid wavelength type: {} (must be array-like)".format(edl.string_prefix,type(wavelength)))

		self.plot_step = self.wavelength[2] - self.wavelength[1]
		#self.__dict__ = dict(np.unique(np.concatenate(([arg for arg in kwargs if arg in edl.keys],[init for init in initial_values])))) # pick up any optionals
		self.caller = 'wavelength'
		self.change()


	''' tier 0 '''

	def __call__(self):
		if (self.initialize_singletons == True):
			self.caller = 'wavelength'
			self.initialize_singletons = False
		else:
			self.caller = cb_obj
		self.change()


	def change(self): # caller=cb_obj.name,tabs=tabs
		caller = self.caller

		# direct the changed values
		if (tabs.active in [_opt for _opt in range(3)]):
			if caller in edl.signal_keys:
				self.recalculate_signal(caller)
			if caller in edl.noise_keys:
				self.recalculate_noise(caller)
			if caller in edl.error_keys:
				self.recalculate_error(caller)
			else:
				self.recalculate_signal(caller)
				self.recalculate_noise(caller)
				self.recalculate_error(caller)
		
			if (tabs.active == 0): # snr
				self.recalculate_snr(caller)
				plot_y1 = self.snr_blue
				plot_y2 = self.snr_red
			elif (tabs.active == 1): # obs spec no noise
				if self.noise:
					plot_y1 = np.add(self.signal_blue,self.error_blue)
					plot_y2 = np.add(self.signal_red,self.error_blue)
				else:
					plot_y1 = self.signal_blue
					plot_y2 = self.signal_red
			elif (tabs.active == 2): # obs sky background
				plot_y1 = self.error_blue
				plot_y2 = self.error_red
			else:
				pass # pass to other options

		else:
			if (tabs.active == 3): # dichroic throughput
				plot_y1 = self.dichro_blue
				plot_y2 = self.dichro_red
			elif (tabs.active == 4): # grating throughput
				plot_y1 = self.grating_blue
				plot_y2 = self.grating_red
			elif (tabs.active == 5): # CCD QE
				plot_y1 = self.ccd_blue
				plot_y2 = self.ccd_red
			elif (tabs.active == 6): # atmospheric extinction
				plot_y1 = self.extinction
				plot_y2 = self.extinction
			else:
				raise ValueError("{} Invalid active tab: {}".format(edl.string_prefix,tabs.active))

			if (cb_obj.name == 'channel'): # change to actual instance's name attribute, not string
				if 0 in cb_obj.active:
					red_active = True
				if 1 in cb_obj.active:
					blue_active = True
				else:
					red_active = True
					blue_active = True
			if red_active and blue_active:
				self.channel = 'both'
			elif blue_active and not red_active:
				self.channel = 'blue'
			elif red_active and not blue_active:
				self.channel = 'red'

		return plot_y1,plot_y2


	def recalculate_snr(self,caller):
		if caller in edl.readnoise_keys:
			self.recalculate_readnoise(caller)
		if caller in edl.signal_keys:
			self.recalculate_signal(caller)
		if caller in edl.noise_keys:
			self.recalculate_noise(caller)
		else:
			self.recalculate_readnoise(caller)
			self.recalculate_signal(caller)
			self.recalculate_noise(caller)

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.snr_blue = np.divide(self.signal_blue,np.sqrt(self.signal_blue + self.noise_blue + np.square(self.readnoise)))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.snr_red = np.divide(self.signal_red,np.sqrt(self.signal_red + self.noise_red + np.square(self.readnoise)))


	def recalculate_error(self,caller):
		if caller in edl.signal_keys:
			self.recalculate_signal(caller)
		if caller in edl.noise_keys:
			self.recalculate_noise(caller)
		if caller in edl.readnoise_keys:
			self.recalculate_readnoise(caller)
		if caller in edl.snr_keys:
			self.recalculate_snr(caller)
		else:
			self.recalculate_signal(caller)
			self.recalculate_noise(caller)
			self.recalculate_readnoise(caller)
			self.recalculate_snr(caller)

		if (self.channel == 'blue') or (self.channel == 'both'):
			sigma_blue = np.sqrt(self.signal_blue + self.noise_blue + np.square(self.readnoise))
			self.error_blue = np.random.normal(loc=0, scale=sigma_blue,size=len(self.snr_blue))
		if (self.channel == 'red') or (self.channel == 'both'):
			sigma_red = np.sqrt(self.signal_red + self.noise_red + np.square(self.readnoise))
			self.error_red = np.random.normal(loc=0, scale=sigma_red,size=len(self.snr_red))


	''' tier 1 '''

	def recalculate_signal(self,caller):
		if caller in edl.signal_keys:
			self.recalculate_signal(caller)
		elif caller in edl.percent_keys:
			self.recalculate_seeing(caller)
		elif caller in edl.total_eff_keys:
			self.recalculate_efficiency(caller)
		else:
			self.recalculate_counts(caller)
			self.recalculate_seeing(caller)
			self.recalculate_efficiency(caller)

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.signal_blue = np.multiply((self.counts * self.percent), self.total_eff_blue)
		if (self.channel == 'red') or (self.channel == 'both'):
			self.signal_red = np.multiply((self.counts * self.percent), self.total_eff_red)

	def recalculate_noise(self,caller):
		if caller in edl.counts_noise_keys:
			self.recalculate_counts_noise(caller)
		else:
			self.recalculate_counts_noise(caller)
			self.recalculalte_efficiency_noise(caller)

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.noise_blue = np.multiply(self.counts_noise,self.total_eff_noise_blue)
		if (self.channel == 'red') or (self.channel == 'both'):
			self.noise_red = np.multiply(self.counts_noise,self.total_eff_noise_red)


	''' tier 2 '''

	def recalculate_counts(self,caller):
		if caller in edl.counts_keys:
			self.recalculate_flux(caller)
			self.change_telescope_mode(caller)
		else:
			self.recalculate_flux(caller)
			self.change_telescope_mode(caller)
			self.change_plot_step(caller)

		self.power = self.flux_y * self.area * self.exposure_time * self.plot_step
		self.counts = np.divide(np.divide(self.power,np.divide((const.h.value * const.c.value),self.wavelength)),1e10)


	def recalculate_counts_noise(self,caller):
		if (caller == 'moon_days'):
			self.recalculate_sky_flux(caller)
		if (caller == 'telescope_mode'):
			self.change_telescope_mode(caller)
		if caller in edl.extension_keys:
			self.recalculate_seeing(caller)
		else:
			self.recalculate_sky_flux(caller)
			self.recalculate_seeing(caller)
			self.change_telescope_mode(caller)
			self.change_plot_step(caller)

		self.counts_noise = np.multiply(np.multiply(self.sky_flux,self.extension),(self.area*self.exposure_time*self.plot_step))


	def recalculate_percent(self,caller):
		self.recalculate_seeing(caller)
		# this is a forwarding function


	def recalculate_efficiency(self,caller):
		if caller in edl.grating_keys:
			self.recalculate_grating(caller)
		else:
			self.recalculate_atmospheric_extinction(caller)
			self.recalculate_grating(caller)
			self.recalculate_dichroic(caller)
			self.recalculate_ccd(caller)
			self.recalculate_atmospheric_extinction(caller)
			self.recalculate_mirror(caller)

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.total_eff_blue = np.multiply(np.multiply(self.dichro_blue,self.grating_blue),np.multiply((self.ccd_blue * (edl.coating_eff_blue * self.extinction)),np.square(self.mirror)))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.total_eff_red = np.multiply(np.multiply(self.dichro_red,self.grating_red),np.multiply((self.ccd_red * (edl.coating_eff_red * self.extinction)),np.square(self.mirror)))


	def recalculalte_efficiency_noise(self,caller):
		if caller in edl.grating_keys:
			self.recalculate_grating(caller)
		else:
			self.recalculate_dichroic(caller)
			self.recalculate_grating(caller)
			self.recalculate_ccd(caller)
			self.recalculate_mirror(caller)

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.total_eff_noise_red = np.multiply(np.multiply(self.dichro_blue,self.grating_blue),(self.ccd_blue * np.square(self.mirror) * edl.coating_eff_blue))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.total_eff_noise_blue = np.multiply(np.multiply(self.dichro_red,self.grating_red),(self.ccd_red * np.square(self.mirror) * edl.coating_eff_red))


	def recalculate_flux(self,caller):
		if caller in edl.grating_opt_keys:
			self.change_grating_opt(caller)
		if caller in edl.filter_keys:
			self.change_filter(caller)
		if (caller == 'object_type'):
			self.change_object_type(caller)
		else: # callers prolly not necessary but good formality
			self.recalculate_atmospheric_extinction(caller)
			self.change_grating_opt(caller)
			self.change_filter(caller)
			self.change_moon_days(caller)
			self.change_object_type(caller)
			self.change_plot_step(caller)

		# heal identicalities
		self.lambda_A[0] = self.lambda_A[0] + self.plot_step
		self.lambda_A[-1] = self.lambda_A[-1] - self.plot_step

		ftrans = interpolate.interp1d(self.selected_filter[0],self.selected_filter[1], kind='cubic')
		self.trans = ftrans(self.lambda_A)

		self._extinction = self.spectres(self.lambda_A,dh.atmo_ext_x,dh.atmo_ext_y)

		self.flux = self.flux_A * 1e10
		self._lambda = self.lambda_A / 1e10

		# valaidate flux
		num_zeros = 0
		for lux in self.flux:
			if (lux is None) or (lux is 0):
				num_zeros += 1
				lux = 0

		if (num_zeros >= (self.flux.shape[0]/5)):
			if (num_zeros == self.flux.shape[0]):
				print('No flux in this bandpass!')
				output_flux = np.zeros(self.wavelength.shape[0])
			else:
				percent_zeros = (num_zeros / self.flux.shape[0]) * 100
				print('{}% of this bandpass has zero flux'.format(percent_zeros))

		self.change_mag_sys_opt(caller)

		del_mag = self.magnitude - self.mag_model
		output_flux = np.multiply(self.object_y,10 ** np.negative(del_mag/2.5))
		old_res = self.object_x[2] - self.object_x[1]
		if (old_res < self.plot_step):
			self.flux_y = self.spectres(self.wavelength,self.object_x,(output_flux*1e-03)) # ergs s-1 cm-2 A-1 to J s-1 m-2 A-1
		else:
			self.flux_y = self.spectres(self.wavelength,self.object_x,(output_flux*1e-03))


	''' tier 3 '''


	def change_grating_opt(self,caller):
		_active = ''
		if self.grating_opt in edl.grating_opt_keys[:2]: # low resolution... or is this backards?
			self.delta_lambda = edl.dld[0] * self.slit_size / 0.7
			_active = 'low-resolution'
		elif self.grating_opt in edl.grating_opt_keys[3:]:
			self.delta_lambda = edl.dld[1] * self.slit_size / 0.7
			_active = 'high-resolution'
		else:
			raise ValueError("{} Invalid grating_opt: {}".format(edl.string_prefix,self.grating_opt))
		if not self.sss:
			print("{} Grating changed to {}".format(edl.string_prefix,_active))


	def change_mag_sys_opt(self,caller):
		try:
			if isinstance(self.old_mso,str):		
				#print('_lambda: {}\ntrans: {}\n_extinction: {}\nflux: {}'.format(type(self._lambda),type(self.trans),type(self._extinction),type(self.flux))) # for debugging
				if (self.mag_sys_opt == 'vega'):
					flux_vega = self.spectres(self.wavelength,dh.vega[0],dh.vega[1]) * 1e10 # fixed... I hope?
					self.mag_model = -2.5 * np.log10(np.divide(math.fsum(self.flux * self._extinction * self._lambda * self.trans),math.fsum(flux_vega * self.trans * self._lambda * self._extinction))) + 0.03
				elif (self.mag_sys_opt == 'ab'):
					sq = np.square(self._lambda)
					self.mag_model = -48.6 - 2.5 * np.log10(math.fsum(self.flux * self.trans * self._extinction * self._lambda) / math.fsum(self.trans * self._lambda * self._extinction * (const.c.value/sq)))
				else:
					raise ValueError("{} Invalid magnitude system option (mag_sys_opt): {}".format(edl.string_prefix,self.mag_sys_opt))
				if not self.sss:
					print("{} Magnitude system changed to {}".format(edl.string_prefix,self.mag_sys_opt.upper()))
				self.old_mso = self.mag_sys_opt
		except:
			pass # no need to update


	def change_object_type(self,caller):
		if self.object_type in edl.stellar_keys:
			index_of = [i for i,name in enumerate(edl.stellar_keys) if self.object_type == name][0]
			self.object_data = dh.starfiles[index_of]
		elif self.object_type in edl.galactic_keys:
			index_of = [i for i,name in enumerate(edl.galactic_keys) if self.object_type == name][0]
			self.object_data = dh.galaxyfiles[index_of]
		else:
			raise ValueError("{} Invalid object type: {}".format(edl.string_prefix,self.object_type))

		self.object_x = self.object_data[0] * (1+self.redshift)
		self.object_y = self.object_data[1]
		self.flux_A = self.spectres(self.lambda_A,self.object_x,self.object_y)
		if not self.sss:
			print("{} Object type changed to {}".format(edl.string_prefix,self.object_type))


	def change_moon_days(self,caller):
		if self.moon_days in edl.moon_days_keys:
			self.sky_background = dh.skyfiles[(int(np.where(np.asarray(edl.moon_days_keys)==self.moon_days)[0]))]
		else:
			raise ValueError('{} Invalid number of days since new moon: {}'.format(edl.string_prefix,self.moon_days))
		self.recalculate_sky_flux(caller)
		if not self.sss:
			print("{} Days since new moon changed to {}".format(edl.string_prefix,self.moon_days))


	def change_telescope_mode(self,caller):
		if self.telescope_mode in edl.telescope_mode_keys[:3]:
			self.area = edl.area[0]
		elif self.telescope_mode in edl.telescope_mode_keys[4:]:
			self.area = edl.area[1]
		else:
			raise ValueError('{} Invalid telescope mode: {}'.format(edl.string_prefix,self.telescope_mode))
		if not self.sss:
			print("{} Telescope mode changed to {} (area: {} m^2)".format(edl.string_prefix,self.telescope_mode,self.area))


	def change_filter(self,caller):
		self.selected_filter = dh.filterfiles[self.filter_index]
		filter_min = min(self.selected_filter[0])
		filter_max = max(self.selected_filter[0])

		if (filter_min > self.wavelength[0]):
			lambda_min = filter_min
		elif (filter_min == self.wavelength[0]):
			filter_min = self.selected_filter[int(np.where(self.selected_filter[0] > self.wavelength[0])[0])]
		else:
			lambda_min = self.wavelength[0]

		if (filter_max < self.wavelength[-1]):
			lambda_max = filter_max
		elif (filter_max == self.wavelength[-1]):
			filter_max = self.selected_filter[int(np.where(self.selected_filter[0] < self.wavelength[-1])[-1])]
		else:
			lambda_max = self.wavelength[-1]

		self.change_plot_step(caller) # i dont remember...
		self.lambda_A = np.arange(lambda_min,lambda_max,self.plot_step)
		if not self.sss:
			_active = edl.filter_files[self.filter_index]
			print("{} Filter changed to {} ({}--{} nm)".format(edl.string_prefix,_active,self.selected_filter[0][0],self.selected_filter[0][-1]))


	def recalculate_seeing(self,caller):
		_sigma = self.seeing / gaussian_sigma_to_fwhm
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		self.percent_u,self.percent_err_u = integrate.quad(funx,(-self.slit_size/2),(self.slit_size/2))
		self.percent_l,self.percent_err_l = integrate.quad(funx,(-self.seeing/2),(self.seeing/2))
		self.percent = self.percent_u * self.percent_l # can use error if you add it later...
		self.extension = self.seeing * self.slit_size
		if not self.sss:
			print("{} Seeing changed to {}".format(edl.string_prefix,self.seeing))


	def recalculate_sky_flux(self,caller):
		sky_x = self.sky_background[0] * 1e4
		sky_y = self.sky_background[1] / 1e4
		old_res = sky_x[2] - sky_x[1]
		_sigma = self.delta_lambda / gaussian_sigma_to_fwhm
		_x = np.arange((-5*_sigma),(5*_sigma),old_res)
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		degrade = funx(_x)/np.trapz(funx(_x))
		sky_y = convolve_fft(sky_y,degrade)
		self.sky_flux = self.spectres(self.wavelength,sky_x,sky_y)


	def recalculate_dichroic(self,caller):
		if (self.channel is 'blue') or (self.channel is 'both'):
			fblue_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y1, kind='cubic')
			self.dichro_blue = fblue_dichro(self.wavelength)
		if (self.channel is 'red') or (self.channel is 'both'):
			fred_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y2, kind='cubic')
			self.dichro_red = fred_dichro(self.wavelength)


	def recalculate_grating(self,caller):
		if (self.channel is 'blue') or (self.channel is 'both'):
			self.grating_blue = self.spectres(self.wavelength,(dh.grating1[0]*10),dh.grating1[1])
		if (self.channel is 'red') or (self.channel is 'both'):
			self.grating_red = self.spectres(self.wavelength,(dh.grating2[0]*10),dh.grating2[1])


	def recalculate_ccd(self,caller):
		if (self.channel is 'blue') or (self.channel is 'both'):
			fblue_ccd = interpolate.interp1d((dh.ccd1[0]*10),dh.ccd1[1], kind='cubic')
			self.ccd_blue = fblue_ccd(self.wavelength)
		if (self.channel is 'red') or (self.channel is 'both'):
			fred_ccd = interpolate.interp1d((dh.ccd2[0]*10),dh.ccd2[1], kind='cubic')
			self.ccd_red = fred_ccd(self.wavelength)

	def recalculate_mirror(self,caller):
		fmirror = interpolate.interp1d(dh.mirror_file_x,dh.mirror_file_y, kind='cubic')
		self.mirror = fmirror(self.wavelength)

	def recalculate_readnoise(self,caller): # probably not implemented in all necessary places yet...
		spectral_resolution = math.ceil((self.slit_size/(0.7/12))/2)*2 # px (ceil()/2)*2 to round up to next even integer
		spatial_resolution = math.ceil((self.seeing/(0.7/12))/2)*2 # probably have to find another method, because pscript is depriciating :(
		extent = self.seeing * self.slit_size
		npix = extent/(0.7/12)**2
		
		rn = edl.rn_default

		if (self.bin_size > 0) and (self.bin_size < 5):
			self.readnoise = math.ceil(rn * spectral_resolution * spatial_resolution / (self.bin_size**2))
			if not self.sss:
				print('{} Pixel binning: ({}x{})'.format(edl.string_prefix,self.bin_size,self.bin_size))
				print('{} Extent: {} arcsec^2\n{} num pixels/resel: {} px\n{} spectral resolution: {} px\n{} spatial resolution: {} px'.format(edl.string_prefix,extent,edl.string_prefix,int(math.ceil(npix)),edl.string_prefix,spectral_resolution,edl.string_prefix,spatial_resolution))
		else:
			raise ValueError('{} Invalid pixel binning option: {}'.format(edl.string_prefix,self.bin_size))


	def recalculate_atmospheric_extinction(self,caller):
		self.extinction = self.spectres(self.wavelength,dh.atmo_ext_x,dh.atmo_ext_y)


	def change_plot_step(self,caller):
		if (self.wavelength[0] >= dfs.default_limits_wavelength[0]) and (self.wavelength[-1] <= dfs.default_limits_wavelength[1]):
			self.plot_step = self.wavelength[2] - self.wavelength[1]
		else: # technically shouldn't happen...
			raise ValueError('{} Invalid wavelength extrema ({}--{}). Must be within {}--{}'.format(edl.string_prefix,wavelength[0],wavelength[-1],dfs.default_limits_wavelength[0],dfs.default_limits_wavelength[1]))


	def num_router(self,_obj):
		try:
			if isinstance(_obj,int):
				return _obj
			elif isinstance(_obj,float):
				return _obj
			else:
				raise ValueError("{} Invalid value: {}".format(edl.string_prefix,_obj))
		except:
			raise TypeError("{} Invalid type: {}".format(edl.string_prefix,type(_obj)))



	def type_router(self,_obj,_keys):
		try:
			if isinstance(_obj,str):
				return int(np.where(np.asarray(_keys)==_obj.lower())[0])
			else:
				return self.num_router(_obj)
		except:
			raise TypeError("{} Invalid type: {}".format(edl.string_prefix,type(_obj)))


	# adaptation of SpectRes, from [Adam Carnall](https://github.com/ACCarnall/SpectRes)
	def spectres(self,new_spec_wavs, old_spec_wavs, spec_fluxes):
	    # Arrays of left-hand sides and widths for the old and new bins
	    spec_lhs = np.zeros(old_spec_wavs.shape[0])
	    spec_widths = np.zeros(old_spec_wavs.shape[0])
	    spec_lhs = np.zeros(old_spec_wavs.shape[0])
	    spec_lhs[0] = old_spec_wavs[0]
	    spec_lhs[0] -= (old_spec_wavs[1] - old_spec_wavs[0])/2
	    spec_widths[-1] = (old_spec_wavs[-1] - old_spec_wavs[-2])
	    spec_lhs[1:] = (old_spec_wavs[1:] + old_spec_wavs[:-1])/2
	    spec_widths[:-1] = spec_lhs[1:] - spec_lhs[:-1]

	    filter_lhs = np.zeros(new_spec_wavs.shape[0]+1)
	    filter_widths = np.zeros(new_spec_wavs.shape[0])
	    filter_lhs[0] = new_spec_wavs[0]
	    filter_lhs[0] -= (new_spec_wavs[1] - new_spec_wavs[0])/2
	    filter_widths[-1] = (new_spec_wavs[-1] - new_spec_wavs[-2])
	    filter_lhs[-1] = new_spec_wavs[-1]
	    filter_lhs[-1] += (new_spec_wavs[-1] - new_spec_wavs[-2])/2
	    filter_lhs[1:-1] = (new_spec_wavs[1:] + new_spec_wavs[:-1])/2
	    filter_widths[:-1] = filter_lhs[1:-1] - filter_lhs[:-2]

	    if filter_lhs[0] < spec_lhs[0] or filter_lhs[-1] > spec_lhs[-1]:
	        raise ValueError("[ SpectRes ] : The new wavelengths specified must fall"
	                         "within the range of the old wavelength values.")

	    # Generate output arrays to be populated
	    res_fluxes = np.zeros(spec_fluxes[0].shape + new_spec_wavs.shape)

	    start = 0
	    stop = 0

	    # Calculate new flux and uncertainty values, loop over new bins
	    for j in range(new_spec_wavs.shape[0]):

	        # Find first old bin which is partially covered by the new bin
	        while spec_lhs[start+1] <= filter_lhs[j]:
	            start += 1

	        # Find last old bin which is partially covered by the new bin
	        while spec_lhs[stop+1] < filter_lhs[j+1]:
	            stop += 1

	        # If new bin is fully within one old bin these are the same
	        if stop == start:
	            res_fluxes[j] = spec_fluxes[start]

	        # Otherwise multiply the first and last old bin widths by P_ij
	        else:
	            start_factor = ((spec_lhs[start+1] - filter_lhs[j])
	                            / (spec_lhs[start+1] - spec_lhs[start]))

	            end_factor = ((filter_lhs[j+1] - spec_lhs[stop])
	                          / (spec_lhs[stop+1] - spec_lhs[stop]))

	            spec_widths[start] *= start_factor
	            spec_widths[stop] *= end_factor

	            # Populate res_fluxes spectrum and uncertainty arrays
	            f_widths = spec_widths[start:stop+1]*spec_fluxes[start:stop+1]
	            res_fluxes[j] = np.sum(f_widths, axis=-1)
	            res_fluxes[j] /= np.sum(spec_widths[start:stop+1])

	            # Put back the old bin widths to their initial values for later use
	            spec_widths[start] /= start_factor
	            spec_widths[stop] /= end_factor

	        return res_fluxes
        # end of simulate.spectres

# end of class `simulate`