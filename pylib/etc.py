#!/usr/bin/python3

''' GMACS Exposure Time Calculator '''

import math
import numpy as np

from spectres import spectres
from scipy import interpolate, integrate
from astropy import constants # `const` not permitted in JS
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_sigma_to_fwhm

import datahandler as dh
import tools as etools
import keys as etkeys
import paths as etpaths
import defaults as dfs


class calculate:
	
	def __init__(self,**kwargs):
		# validate request
		self.__dict__ = kwargs
		self.wavelength = np.asarray(self.wavelength)
		self.mag_sys_opt = self.mag_sys # fix this later
		if self.telescope_mode == 'full':
			self.num_mirrors = 7
		else:
			self.num_mirrors = 4
		#self.__dict__,self.wavelength = etools.validate_args(kwargs)
		self.update()

		
	def update(self):
		if (self.mode == 0):
			plot_y = self.snr()
		elif (self.mode == 1):
			plot_y = self.os()
		elif (self.mode == 2):
			plot_y = self.sky()
		elif (self.mode == 3):
			plot_y = self.dichro()
		elif (self.mode == 4):
			plot_y = self.gt()
		elif (self.mode == 5):
			plot_y = self.ccd()
		elif (self.mode == 6):
			plot_y = self.atmo_ext()
		else:
			raise ValueError("{} Invalid calculator mode: ({})".format(dfs.string_prefix,self.mode))
		if (self.mode == 6): # exception for atmo ext
			values_y_blue = plot_y[0]
			values_y_red = plot_y[1]
		else:
			#values_x,values_y_blue,values_y_red = etools.validate_points(self.wavelength,plot_y[0])
			values_y_blue = plot_y[:int(len(plot_y)/2)]
			values_y_red = plot_y[int(len(plot_y)/2):]

		self.values = [self.wavelength,values_y_blue,values_y_red]
		

	# mode functions
	def snr(self):
		plot_y = np.empty(0)
		self.calc_snr()
		if (self.channel == 'blue') or (self.channel == 'both'):
			plot_y = np.append(plot_y,self.snr_blue)
		if (self.channel == 'red') or (self.channel == 'both'):
			plot_y = np.append(plot_y,self.snr_red)
		else:
			raise ValueError("{} Invalid channel selection: ({})".format(dfs.string_prefix,self.channel))
		return plot_y

	def os(self):
		plot_y = np.empty(0)
		calc_error()
		if self.noise:
			if (self.channel == 'blue') or (self.channel == 'both'):
				plot_y = np.append(plot_y,np.add(self.signal_blue,self.error_blue))
			if (self.channel == 'red') or (self.channel == 'both'):
				plot_y = np.append(plot_y,np.add(self.signal_red,self.error_red))
			else:
				raise ValueError("{} Invalid channel selection: ({})".format(dfs.string_prefix,self.channel))
		else:
			if (self.channel == 'blue') or (self.channel == 'both'):
				plot_y1 = self.signal_blue
			if (self.channel == 'red') or (self.channel == 'both'):
				plot_y = np.append(plot_y,self.signal_red)
			else:
				raise ValueError("{} Invalid channel selection: ({})".format(dfs.string_prefix,self.channel))
		return plot_y

	def sky(self):
		plot_y = np.empty(0)
		self.calc_err()
		if (self.channel == 'blue') or (self.channel == 'both'):
			plot_y = np.append(plot_y,self.error_blue)
		if (self.channel == 'red') or (self.channel == 'both'):
			plot_y = np.append(plot_y,self.error_red)
		else:
			raise ValueError("{} Invalid channel selection: ({})".format(dfs.string_prefix,self.channel))
		return plot_y

	def dichro(self):
		plot_y = np.empty(0)
		self.recalculate_dichroic()
		if (self.channel == 'blue') or (self.channel == 'both'):
			plot_y = np.append(plot_y,self.dichro_blue)
		if (self.channel == 'red') or (self.channel == 'both'):
			plot_y = np.append(plot_y,self.dichro_red)
		else:
			raise ValueError("{} Invalid channel selection: ({})".format(dfs.string_prefix,self.channel))
		return plot_y

	def gt(self):
		plot_y = np.empty(0)
		self.recalculate_grating()
		if (self.channel == 'blue') or (self.channel == 'both'):
			plot_y = np.append(plot_y,self.grating_blue)
		if (self.channel == 'red') or (self.channel == 'both'):
			plot_y = np.append(plot_y,self.grating_red)
		else:
			raise ValueError("{} Invalid channel selection: ({})".format(dfs.string_prefix,self.channel))
		return plot_y

	def ccd(self):
		plot_y = np.empty(0)
		self.recalculate_ccd()
		if (self.channel == 'blue') or (self.channel == 'both'):
			plot_y = np.append(plot_y,self.ccd_blue)
		if (self.channel == 'red') or (self.channel == 'both'):
			plot_y = np.append(plot_y,self.ccd_red)
		else:
			raise ValueError("{} Invalid channel selection: ({})".format(dfs.string_prefix,self.channel))
		return plot_y

	def atmo_ext(self):
		self.recalculate_atmospheric_extinction()
		return self.extinction

	def calc_snr(self):
		self.recalculate_readnoise()
		self.recalculate_signal()
		self.recalculate_noise()
		if (self.channel == 'blue') or (self.channel == 'both'):
			self.snr_blue = np.divide(self.signal_blue,np.sqrt(self.signal_blue + self.noise_blue + np.square(self.read_noise)))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.snr_red = np.divide(self.signal_red,np.sqrt(self.signal_red + self.noise_red + np.square(self.read_noise)))

	def calc_error(self):
		self.recalculate_readnoise()
		self.recalculate_signal()
		self.recalculate_noise()
		if (self.self.channel == 'blue') or (self.channel == 'both'):
			sigma_blue = np.sqrt(self.signal_blue + self.noise_blue + np.square(self.read_noise))
			self.error_blue = np.random.normal(loc=0, scale=sigma_blue,size=len(self.snr_blue))
		if (self.channel == 'red') or (self.channel == 'both'):
			sigma_red = np.sqrt(self.signal_red + self.noise_red + np.square(self.read_noise))
			self.error_red = np.random.normal(loc=0, scale=sigma_red,size=len(self.snr_red))

	def recalculate_signal(self):
		self.recalculate_counts()
		self.recalculate_seeing()
		self.recalculate_efficiency()
		if (self.channel == 'blue') or (self.channel == 'both'):
			self.signal_blue = np.multiply((self.counts * self.percent), self.total_eff_blue)
		if (self.channel == 'red') or (self.channel == 'both'):
			self.signal_red = np.multiply((self.counts * self.percent), self.total_eff_red)

	def recalculate_noise(self):
		self.recalculate_counts_noise()
		self.recalculalte_efficiency_noise()
		if (self.channel == 'blue') or (self.channel == 'both'):
			self.noise_blue = np.multiply(self.counts_noise,self.total_eff_noise_blue)
		if (self.channel == 'red') or (self.channel == 'both'):
			self.noise_red = np.multiply(self.counts_noise,self.total_eff_noise_red)

	def recalculate_counts(self):
		self.recalculate_plot_step()
		self.recalculate_flux()
		self.change_num_mirrors()
		self.power = self.flux_y * self.area * self.exposure_time * self.plot_step
		self.counts = np.divide(np.divide(self.power,np.divide((constants.h.value * constants.c.value),self.wavelength)),1e10)

	def recalculate_counts_noise(self):
		self.recalculate_plot_step()
		self.recalculate_sky_flux()
		self.recalculate_seeing()
		self.change_num_mirrors()
		self.counts_noise = np.multiply(np.multiply(self.sky_flux,self.extension),(self.area*self.exposure_time*self.plot_step))

	def recalculate_percent(self):
		self.recalculate_seeing()
		# this is a forwarding function

	def recalculate_efficiency(self):
		self.recalculate_grating()
		self.recalculate_dichroic()
		self.recalculate_ccd()
		self.recalculate_atmospheric_extinction()
		self.recalculate_mirror_loss()

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.total_eff_blue = np.multiply(np.multiply(self.dichro_blue,self.grating_blue),np.multiply((self.ccd_blue * (dfs.coating_eff_blue * self.extinction)),np.square(self.mirror)))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.total_eff_red = np.multiply(np.multiply(self.dichro_red,self.grating_red),np.multiply((self.ccd_red * (dfs.coating_eff_red * self.extinction)),np.square(self.mirror)))

	def recalculalte_efficiency_noise(self):
		self.recalculate_dichroic()
		self.recalculate_grating()
		self.recalculate_ccd()
		self.recalculate_mirror_loss()
		if (self.channel == 'blue') or (self.channel == 'both'):
			self.total_eff_noise_red = np.multiply(np.multiply(self.dichro_blue,self.grating_blue),(self.ccd_blue * np.square(self.mirror) * dfs.coating_eff_blue))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.total_eff_noise_blue = np.multiply(np.multiply(self.dichro_red,self.grating_red),(self.ccd_red * np.square(self.mirror) * dfs.coating_eff_red))

	def recalculate_flux(self):
		self.recalculate_plot_step()
		self.change_grating_opt()
		self.change_moon_days()
		self.change_filter()
		self.change_object_type()

		# heal identicalities
		self.lambda_A[0] = self.lambda_A[0] + self.plot_step
		self.lambda_A[-1] = self.lambda_A[-1] - self.plot_step

		ftrans = interpolate.interp1d(self.selected_filter[0],self.selected_filter[1], kind='cubic')
		self.trans = ftrans(self.lambda_A)
		self._extinction = spectres(self.lambda_A,dh.atmo_ext_x,dh.atmo_ext_y)
		
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

		self.change_mag_sys_opt()
		del_mag = self.mag - self.mag_model
		output_flux = np.multiply(self.object_y,10 ** np.negative(del_mag/2.5))
		old_res = self.object_x[2] - self.object_x[1]
		if (old_res < self.plot_step):
			self.flux_y = spectres(self.wavelength,self.object_x,(output_flux*1e-03)) # ergs s-1 cm-2 A-1 to J s-1 m-2 A-1
		else:
			self.flux_y = spectres(self.wavelength,self.object_x,(output_flux*1e-03))
		
	def recalculate_plot_step(self):
		if (self.wavelength[0] >= dfs.wavelength_limits[0]) and (self.wavelength[-1] <= dfs.wavelength_limits[1]):
			self.plot_step = self.wavelength[2] - self.wavelength[1]
		else: # technically shouldn't happen...
			raise ValueError('{} Invalid wavelength extrema ({}--{}). Must be within {}--{}'.format(dfs.string_prefix,self.wavelength[0],self.wavelength[-1],dfs.wavelength_limits[0],dfs.wavelength_limits[1]))

	def change_grating_opt(self):
		self.delta_lambda = dfs.dld[self.grating_opt] * self.slit_width / 0.7
		_active = 'high' if (self.grating_opt==1) else 'low'
		if not self.sss:
			print("{} Grating changed to \'{} resolution\'".format(dfs.string_prefix,_active))

	def change_mag_sys_opt(self):
		if (self.mag_sys_opt == 'vega') or (self.mag_sys_opt ==0):
			flux_vega = spectres(self.wavelength,dh.vega[0],dh.vega[1]) * 1e10 # fixed... I hope?
			self.mag_model = -2.5 * np.log10(np.divide(math.fsum(self.flux * self._extinction * self._lambda * self.trans),math.fsum(self.flux_vega * self.trans * self._lambda * self._extinction))) + 0.03
		elif (self.mag_sys_opt == 'ab') or (self.mag_sys_opt ==1):
			self.mag_model = -48.6 - 2.5 * np.log10(math.fsum(self.flux * self.trans * self._extinction * self._lambda) / math.fsum(self.trans * self._lambda * self._extinction * (constants.c.value/np.square(self._lambda))))
		else:
			raise ValueError("{} Invalid magnitude system option (mag_sys_opt): {}".format(dfs.string_prefix,self.mag_sys_opt))
		if not self.sss:
			print("{} Magnitude system changed to {}".format(dfs.string_prefix,np.asarray(etkeys.mag_sys).flatten()[self.mag_sys_opt]))

	def change_object_type(self):
		#self.object_type = dh.galaxyfiles[self.object_type[0]][self.object_type[1]]
		if self.object_class in etkeys.stellar:
			self.object_type = dh.starfiles[int(np.where(np.asarray(etkeys.stellar)==self.object_class)[0])]
		elif self.object_class in etkeys.galactic:
			self.object_type = dh.galaxyfiles[int(np.where(np.asarray(etkeys.stellar)==self.object_class)[0])]
		else:
			raise ValueError("{} Invalid object class: ({})".format(dfs.string_prefix,self.object_class))
		self.object_x = self.object_type[0] * (1 + self.redshift)
		self.object_y = self.object_type[1]
		self.flux_A = spectres(self.lambda_A,self.object_x,self.object_y)
		if not self.sss:
			#print("{} Object type changed to {}".format(dfs.string_prefix,self.object_type))
			pass

	def change_moon_days(self):
		if self.moon_days in range(len(etkeys.moon_days)):
			self.sky_background = dh.skyfiles[(int(np.where(np.asarray(etkeys.moon_days)==self.moon_days)[0]))]
		else:
			raise ValueError('{} Invalid number of days since new moon: {}'.format(dfs.string_prefix,self.moon_days))
		self.recalculate_sky_flux()
		if not self.sss:
			print("{} Days since new moon changed to {}".format(dfs.string_prefix,self.moon_days))

	def change_num_mirrors(self):
		if (self.num_mirrors == 0) or (self.num_mirrors == 4):
			self.area = dfs.area[0]
		elif (self.num_mirrors == 1) or (self.num_mirrors == 7):
			self.area = dfs.area[1]
		else:
			raise ValueError("{} Invalid telescope mode: ({})".format(dfs.string_prefix,self.num_mirrors))
		if not self.sss:
			print("{} Telescope mode changed to: {} m^2".format(dfs.string_prefix,self.area))

	def change_filter(self):
		self.selected_filter = dh.filterfiles[self.filter_opt]
		filter_min = min(self.selected_filter[0])
		filter_max = max(self.selected_filter[0])
		lambda_min = filter_min
		lambda_max = filter_max

		self.recalculate_plot_step()
		self.lambda_A = np.arange(lambda_min,lambda_max,self.plot_step)
		if not self.sss:
			_active = etpaths.filter_files[self.filter_opt]
			print("{} Filter changed to {} ({}--{} nm)".format(dfs.string_prefix,etkeys.filter_opt[0][self.filter_opt],self.selected_filter[0][0],self.selected_filter[0][-1]))

	def recalculate_seeing(self):
		_sigma = self.seeing / gaussian_sigma_to_fwhm
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		self.percent_u,self.percent_err_u = integrate.quad(funx,(-self.slit_width/2),(self.slit_width/2))
		self.percent_l,self.percent_err_l = integrate.quad(funx,(-self.seeing/2),(self.seeing/2))
		self.percent = self.percent_u * self.percent_l # can use error if you add it later...
		self.extension = self.seeing * self.slit_width
		if not self.sss:
			print("{} Seeing changed to {}".format(dfs.string_prefix,self.seeing))

	def recalculate_sky_flux(self):
		sky_x = self.sky_background[0] * 1e4
		sky_y = self.sky_background[1] / 1e4
		old_res = sky_x[2] - sky_x[1]
		_sigma = self.delta_lambda / gaussian_sigma_to_fwhm
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		_x = np.arange((-5*_sigma),(5*_sigma),old_res)
		degrade = funx(_x)/np.trapz(funx(_x))
		sky_y = convolve_fft(sky_y,degrade)
		self.sky_flux = spectres(self.wavelength,sky_x,sky_y)

	def recalculate_dichroic(self):
		if (self.channel is 'blue') or (self.channel is 'both'):
			fblue_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y1, kind='cubic')
			self.dichro_blue = fblue_dichro(self.wavelength)
		if (self.channel is 'red') or (self.channel is 'both'):
			fred_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y2, kind='cubic')
			self.dichro_red = fred_dichro(self.wavelength)

	def recalculate_grating(self):
		if (self.channel is 'blue') or (self.channel is 'both'):
			self.grating_blue = spectres(self.wavelength,(dh.grating1[0]*10),dh.grating1[1])
		if (self.channel is 'red') or (self.channel is 'both'):
			self.grating_red = spectres(self.wavelength,(dh.grating2[0]*10),dh.grating2[1])

	def recalculate_ccd(self):
		if (self.channel is 'blue') or (self.channel is 'both'):
			fblue_ccd = interpolate.interp1d((dh.ccd1[0]*10),dh.ccd1[1], kind='cubic')
			self.ccd_blue = fblue_ccd(self.wavelength)
		if (self.channel is 'red') or (self.channel is 'both'):
			fred_ccd = interpolate.interp1d((dh.ccd2[0]*10),dh.ccd2[1], kind='cubic')
			self.ccd_red = fred_ccd(self.wavelength)

	def recalculate_mirror_loss(self):
		fmirror = interpolate.interp1d(dh.mirror_file_x,dh.mirror_file_y, kind='cubic')
		self.mirror = fmirror(self.wavelength)

	def recalculate_readnoise(self): # probably not implemented in all necessary places yet...
		spectral_resolution = math.ceil((self.slit_width/(0.7/12))/2)*2 # px (ceil()/2)*2 to round up to next even integer
		spatial_resolution = math.ceil((self.seeing/(0.7/12))/2)*2 # probably have to find another method, because pscript is depriciating :(
		extent = self.seeing * self.slit_width
		npix = extent/(0.7/12)**2
		
		rn = dfs.rn_default

		if (self.bin_size > 0) and (self.bin_size < 5):
			self.read_noise = math.ceil(rn * spectral_resolution * spatial_resolution / (self.bin_size**2))
			if not self.sss:
				print('{0} Pixel binning: ({1}x{1})'.format(dfs.string_prefix,self.bin_size))
				print('{0} Extent: {1} arcsec^2\n{0} num pixels/resel: {2} px\n{0} spectral resolution: {3} px\n{0} spatial resolution: {4} px'.format(dfs.string_prefix,extent,int(math.ceil(npix)),spectral_resolution,spatial_resolution))
		else:
			raise ValueError('{} Invalid pixel binning option: {}'.format(dfs.string_prefix,self.bin_size))

	def recalculate_atmospheric_extinction(self):
		self.extinction = spectres(self.wavelength,dh.atmo_ext_x,dh.atmo_ext_y)