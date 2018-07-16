#!/usr/bin/python3
# bork bork bork

import math
import numpy as np

class bork:
	def __init__(self,mag_sys_opt,mag,grating_type,redshift,object_type,filter_opt,seeing,slit_size,area,start_lambda,
		end_lambda,plot_step,exp_time,active_tab,wavelength,colors,sky_files):
		# inits for the wealthy barber (TODO this could probably be more efficient by replacement rather than re-assignment)
		self.mag_sys_opt = mag_sys_opt
		self.mag = mag
		self.grating_type = grating_type
		self.redshift = redshift
		self.object_type = object_type
		self.filter_opt = filter_opt
		self.seeing = seeing
		self.slit_size = slit_size
		self.area = area
		self.start_lambda = start_lambda
		self.end_lambda = end_lambda
		self.plot_step = plot_step
		self.exp_time = exp_time
		self.active_tab = active_tab
		self.wavelength = wavelength
		self.colors = len(colors)
		# specifics
		self.sky_x = sky_files['a']
		self.sky_y = sky_files['b']
		# constants
		self.planck = 6.626 * 10**(-34)						# reduced planck's
		self.celerity = 29979245							# speed of light
		self.hc = self.planck * self.celerity				# conservationism
		self.fwhm = 2.35482 								# full width at half maximum
		setattr(self,'coating_efficiency', 0.8*0.98**14)	# co-efficient of the mirror, plus 14 lens
															# just wanted to see if setattr behaved like that
	def get_delta_lambda(self):
		if (self.grating_type == 0):
			self.grating_opt_val = 3.73
		elif (self.grating_type == 1):
			self.grating_opt_val = 1.40
		else:
			print("[Error] :: Delta Lambda")
		self.delta_lambda = self.grating_opt_val*self.slit_size/0.7

	def get_wavelength(self):
		self.get_delta_lambda()
		_wavelength = np.arange(self.start_lambda,self.end_lambda,self.plot_step)
		self.wavelength = _wavelength[:,np.newaxis]
		return self.wavelength

	def get_sigma(self):
		self.sigma = np.divide(self.seeing,self.fwhm)
		return self.sigma

	def get_power(self):
		self.flux_y = interpolate.PchipInterpolator(self.star_x,self.star_y)(get_wavelength())
		_power = (np.dot(self.flux_y,1e-03).dot(self.area).dot(self.exp_time).dot(self.plot_step))
		self.power = _power[:,np.newaxis]
		return self.power

	def get_counts(self):
		self.counts = np.divide(get_power(),np.divide((self.hc),get_wavelength())) # cuz these are the same
		return self.counts

	# sky land (it's a mario reference)
	def get_sky_counts_noise(self):
		self.old_res = self.sky_x[1] - self.sky_x[0]										# i think i used to have a voice
		self.sky_sigma = np.divide(self.delta_lambda,self.fwhm)								# now i never make a sound
		self.sky_extrema = np.dot(5,self.sky_sigma)											# i just do what i've been told
		_x = np.arange((-self.sky_extrema),(self.sky_extrema),self.old_res)					# i really don't want
		__funx = np.dot(self.sky_sigma,self.sky_sigma).dot(2)
		___funx = np.divide(math.exp())

		self.sky_funx = np.divide()
		self.sky_funx = (1/(self.sky_sigma*math.sqrt(2*math.pi)))*(math.exp((-_x**2)/(2*self.sky_sigma**2)))
		self.sky_degrade = self.sky_funx/np.trapz(self.sky_funx)
		self.sky_sky_y = signal.convolve2d(self.sky_y,self.sky_degrade,'same') 				# convolution to degreade the spectrum

		self.sky_flux = interpolate.PchipInterpolator(self.sky_x,self.sky_sky_y,self.wavelength) # photons/s/m^2/A/arcsec^2
		self.sky_red = np.dot(self.sky_flux,self.extension).dot(self.area).dot(self.exp_time).dot(self.plot_step) # n photons
		self.counts_noise_red = self.sky_red[:,np.newaxis]
		self.counts_noise_blue = self.counts_noise_red
		# return self.sky_red, self.sky_blue, self.counts_noise_red, self.counts_noise_blue

	def get_percent(self):
		self.get_sigma()
		self.funx = (np.divide(1,np.dot(self.sigma,(math.sqrt(2*math.pi))))**((-x**2)/(np.dot(2,self.sigma**2))))
		self.slit_extrema = np.divide(self.slit_size,2)
		self.seeing_extrema = np.divide(self.seeing,2)
		self.integral1,_error = integrate.quad(self.funx,-self.slit_extrema,self.slit_extrema) # error currently unused, figured good to be aware of tho
		self.integral2,_error = integrate.quad(self.funx,-self.seeing_extrema,self.seeing_extrema)
		self.percent = np.dot(self.integral1,self.integral2)

# dichroic transmission
	def get_dichro(self):
		if (self.colors == 0 or self.colors == 2):#### THIS IS WORNG STUFF BELOW@#####
			self.dichro_red = interpolate.PchipInterpolator(_dichro_x,_dichro_y1)(self.wavelength)
		if (self.colors == 1 or self.colors == 2):
			self.dichro_blue = interpolate.PchipInterpolator(_dichro_x,_dichro_y)(self.wavelength)

	def get_signal(self):
		self.get_dichro()
		if (self.colors == 0 or self.colors == 2):
			self.total_efficiency_red = np.dot(self.dichro_red,self.grating_red,self.ccd_red,self.coating_efficiency,self.extinction)
			self.signal_red = np.dot(self.counts,self.percent,self.total_efficiency_red)
		if (self.colors == 1 or self.colors == 2):
			self.total_efficiency_blue = np.dot(self.dichro_blue,self.grating_blue,self.ccd_blue,self.coating_efficiency,self.extinction)
			self.signal_blue = np.dot(self.counts,self.percent,self.total_efficiency_blue)

	def get_noise(self):
		self.get_sky_counts_noise()
		self.get_signal()
		if (self.colors == 0 or self.colors == 2):
			self.total_efficiency_noise_red = np.dot(self.dichro_red,self.grating_red,self.ccd_red,self.coating_efficiency)
			self.noise_red = np.dot(self.counts_noise_red, self.total_efficiency_noise_red)
		if (self.colors == 1 or self.colors == 2):
			self.total_efficiency_noise_blue = np.dot(self.dichro_blue,self.grating_blue,self.ccd_blue,self.coating_efficiency)
			self.noise_blue = np.dot(self.counts_noise_blue, self.total_efficiency_noise_blue)

	def get_snr(self):
		self.get_signal()
		self.get_noise()
		if (self.colors == 0 or self.colors == 2):
			_snr_red = self.signal_red / math.sqrt(self.signal_red + self.noise_red)
			self.snr_red = _snr_red[:,np.newaxis]
		if (self.colors == 1 or self.colors == 2):
			_snr_blue = self.signal_blue / math.sqrt(self.signal_blue + self.noise_blue)
			self.snr_blue = _snr_blue[:,np.newaxis]

	def snr_based_error(self):
		self.get_snr()
		if (self.colors == 0 or self.colors == 2):
			self.snr_based_sigma_red = math.sqrt(self.signal_red+self.noise_red)
			self.snr_based_error_red = np.random.normal(0,self.snr_based_sigma_red)
		if (self.colors == 1 or self.colors == 2):
			self.snr_based_sigma_blue = math.sqrt(self.signal_blue+self.noise_blue)
			self.snr_based_error_blue = np.random.normal(0,self.snr_based_sigma_blue)

	def update(self,widget):
		pass