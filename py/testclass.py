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

	def __init__(self,telescope_mode='first',wavelength=np.arange(3200,10360,edl.dld[0]/3.),exposure_time=3600,object_type='a5v',
		filter_index=3,mag_sys_opt='ab',magnitude=25,redshift=0,seeing=0.5,slit_size=0.5,moon_days=0,grating_opt=0,noise=False,channel='both'):
		object_type = type_router(object_type)
		grating_opt = type_router(grating_opt)
		telescope_mode = type_router(telescope_mode)
		plot_step = wavelength[2] - wavelength[1]


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
		pass


	def recalculate_error(self,caller):
		pass


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


	def recalculate_noise(self,caller):
		if caller in edl.counts_noise_components:
			pass # recalc counts noise stuff
		elif caller in edl.efficiency_noise_components:
			pass # recalc total efficiency noise stuff


	''' tier 2 '''

	def recalculate_counts(self,caller):
		flux = self.flux_y
		self.power = flux * area * exp_time * plot_step
		self.counts = np.divide(np.divide(self.power,np.divide((const.h.value * const.c.value),self.wavelength)),1e10)


	def recalculate_counts_noise(self,caller):
		recalculate_sky_flux(caller)
		recalculate_extension(caller)
		self.counts_noise = np.multiply(np.multiply(self.sky_flux,self.extension),(self.area*self.exp_time*self.plot_step))


	def recalculate_percent(self,caller):
		pass # recalc percent stuff


	def recalculate_efficiency(self,caller):
		pass # recalc efficiency stuff

	def recalculalte_efficiency_noise(self,caller):
		pass # recalc efficienccy noise stuff


	def recalculate_flux(self,caller):
		# heal identicalities
		lambda_A[0] = lambda_A[0] + plot_step
		lambda_A[-1] = lambda_A[-1] - plot_step

		ftrans = interpolate.interp1d(selected_filter[0],selected_filter[1], kind='cubic')
		trans = ftrans(lambda_A) #spectres(lambda_A,selected_filter[0],selected_filter[1])

		extinction = spectres(lambda_A,atmo_ext_x,atmo_ext_y)

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
			mag_model = -2.5 * np.log10(np.divide(math.fsum(flux * extinction * _lambda * trans),math.fsum(flux_vega * trans * _lambda * extinction))) + 0.03
		elif (mag_sys_opt == 'ab'):
			#mag_model = -48.60 - 2.5 * np.log10(np.divide(math.fsum(np.multiply(flux,trans).multiply(extinction).multiply(_lambda)),math.fsum(trans.multiply(_lambda) * extinction[1]).multiply(const.c.value/(np.square(_lambda)))))
			mag_model = -48.6 - 2.5 * np.log10(math.fsum(flux * trans * extinction *_lambda) / math.fsum(trans * _lambda * extinction * (const.c.value/np.square(_lambda))))
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

