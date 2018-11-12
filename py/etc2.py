import math
import numpy as np

import spectres as spectres
from scipy import interpolate, integrate
from astropy import constants # `const` not permitted in JS
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_sigma_to_fwhm

import tools as etools
import keys as etkeys
import paths as etpaths
import defaults as dfs


class calculate:
	
	def __init__(self,**kwargs):
		# validate request
		for _key,_value in kwargs.items():
			validated_args,validated_values = np.empty(0),np.empty(0)
			if _key in np.asarray(etkeys.arguments).flatten():
				validated_args = np.append(validated_args,etools.arg_registrar(_key))
			else:
				raise ValueError("{} Invalid argument: ({})".format(dfs.string_prefix,_key))
			
			if _value in np.asarray(etkeys.keychain).flatten():
				validated_values = np.append(validated_values,etools.value_registrar(_value))
			else:
				raise ValueError("{} Invalid value: ({})".format(dfs.string_prefix,_value))

		if (len(validated_args) != len(validated_values)):
			raise ValueError("{} Argument / value mis-match!".format(dfs.string_prefix))
		else:
			self.__dict__ = zip(validated_args,validated_values)

		# run initial calculation
		if (self.mode == 0):
			snr()

		elif (self.mode == 1):
			os()

		elif (self.mode == 2):
			# sky background
		elif (self.mode == 3):
			# dichroic
		else:
			# etc

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

	def atmo_extinction(self):
		self.recalculate_atmospheric_extinction()
		return [self.extinction]

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
		if (channel == 'red') or (self.channel == 'both'):
			self.noise_red = np.multiply(self.counts_noise,self.total_eff_noise_red)

	def recalculate_counts(self):
		self.recalculate_flux()
		self.change_telescope_mode()
		self.change_plot_step()
		self.power = self.flux_y * self.area * self.exposure_time * self.plot_step
		self.counts = np.divide(np.divide(self.power,np.divide((constants.h.value * constants.c.value),self.wavelength)),1e10)

	def recalculate_counts_noise(self):
		self.recalculate_sky_flux()
		self.recalculate_extension()
		self.change_telescope_mode()
		self.change_plot_step()
		self.counts_noise = np.multiply(np.multiply(self.sky_flux,self.extension),(self.area*self.exposure_time*self.plot_step))

	def recalculate_percent(self):
		self.recalculate_seeing()
		# this is a forwarding function

	def recalculate_efficiency(self):
		self.recalculate_grating()
		self.recalculate_dichroic()
		self.recalculate_ccd()
		self.recalculate_atmospheric_extinction()
		self.recalculate_mirror()

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.total_eff_blue = np.multiply(np.multiply(self.dichro_blue,self.grating_blue),np.multiply((self.ccd_blue * (edl.coating_eff_blue * self.extinction)),np.square(mirror)))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.total_eff_red = np.multiply(np.multiply(self.dichro_red,self.grating_red),np.multiply((self.ccd_red * (edl.coating_eff_red * self.extinction)),np.square(mirror)))

	def recalculalte_efficiency_noise(self):
		self.recalculate_dichroic()
		self.recalculate_grating()
		self.recalculate_ccd()
		self.recalculate_mirror()
		if (self.channel == 'blue') or (self.channel == 'both'):
			self.total_eff_noise_red = np.multiply(np.multiply(self.dichro_blue,self.grating_blue),(self.ccd_blue * np.square(self.mirror) * dfs.coating_eff_blue))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.total_eff_noise_blue = np.multiply(np.multiply(self.dichro_red,self.grating_red),(self.ccd_red * np.square(self.mirror) * dfs.coating_eff_red))

	def recalculate_flux(self):
		self.change_object_type()
		self.change_grating_opt()
		self.change_filter()
		self.change_moon_days()
		self.change_plot_step()

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
		for lux in flux:
			if (lux is None) or (lux is 0):
				num_zeros += 1
				lux = 0

		if (num_zeros >= (flux.shape[0]/5)):
			if (num_zeros == flux.shape[0]):
				print('No flux in this bandpass!')
				output_flux = np.zeros(self.wavelength.shape[0])
			else:
				percent_zeros = (num_zeros / flux.shape[0]) * 100
				print('{}% of this bandpass has zero flux'.format(percent_zeros))

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
		_active = ['high' if (self.grating_opt==1) else 'low']
		if not self.sss:
			print("{} Grating changed to {}".format(dfs.string_prefix,_active))

	def change_mag_sys_opt(self):
		if (self.mag_sys_opt == 'vega'):
			flux_vega = spectres(self.wavelength,dh.vega[0],dh.vega[1]) * 1e10 # fixed... I hope?
			self.mag_model = -2.5 * np.log10(np.divide(math.fsum(self.flux * self._extinction * self._lambda * self.trans),math.fsum(self.flux_vega * self.trans * self._lambda * self._extinction))) + 0.03
		elif (self.mag_sys_opt == 'ab'):
			self.mag_model = -48.6 - 2.5 * np.log10(math.fsum(self.flux * self.trans * self._extinction * self._lambda) / math.fsum(self.trans * self._lambda * self._extinction * (constants.c.value/np.square(self._lambda))))
		else:
			raise ValueError("{} Invalid magnitude system option (mag_sys_opt): {}".format(dfs.string_prefix,self.mag_sys_opt))
		if not self.sss:
			print("{} Magnitude system changed to {}".format(dfs.string_prefix,self.mag_sys_opt.upper()))

	def change_object_type(self):
		if self.object_type in etkeys.stellar:
			index_of = [i for i,name in enumerate(etkeys.stellar) if self.object_type in name][0]
			self.object_type = dh.starfiles[index_of]
		elif self.object_type in etkeys.galactic:
			index_of = [i for i,name in enumerate(etkeys.galactic) if self.object_type in name][0]
			self.object_type = dh.galaxyfiles[index_of]
		else:
			raise ValueError("{} Invalid object type: {}".format(dfs.string_prefix,self.object_type))
		self.object_x = self.object_type[0] * (1 + self.redshift)
		self.object_y = self.object_type[1]
		self.flux_A = spectres(self.lambda_A,self.object_x,self.object_y)
		if not self.sss:
			print("{} Object type changed to {}".format(dfs.string_prefix,self.object_type))

	def change_moon_days(self):
		if self.moon_days in etkeys.moon_days:
			self.sky_background = dh.skyfiles[(int(np.where(np.asarray(etkeys.moon_days)==self.moon_days)[0]))]
		else:
			raise ValueError('{} Invalid number of days since new moon: {}'.format(dfs.string_prefix,self.moon_days))
		self.recalculate_sky_flux(caller)
		if not self.sss:
			print("{} Days since new moon changed to {}".format(dfs.string_prefix,self.moon_days))

	def change_num_mirrors(self):
		if (self.num_mirrors == 0):
			self.area = dfs.area[0]
		elif (self.num_mirrors == 1):
			self.area = dfs.area[1]
		else:
			raise ValueError("{} Invalid telescope mode: ({})".format(dfs.string_prefix,self.num_mirrors))
		if not self.sss:
			print("{} Telescope mode changed to: {} m^2".format(dfs.string_prefix,self.area))

	def change_filter_opt(self):
		self.selected_filter = dh.filterfiles[self.filter_opt]
		filter_min = min(self.selected_filter[0])
		filter_max = max(self.selected_filter[0])
		lambda_min = filter_min
		lambda_max = filter_max

		self.recalculate_plot_step()
		self.lambda_A = np.arange(lambda_min,lambda_max,self.plot_step)
		if not self.sss:
			_active = etpaths.filter_files[self.filter_index]
			print("{} Filter changed to {} ({}--{} nm)".format(dfs.string_prefix,_active,self.selected_filter[0],self.selected_filter[-1]))

	def recalculate_seeing(self):
		_sigma = self.seeing / gaussian_sigma_to_fwhm
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		self.percent_u,self.percent_err_u = integrate.quad(funx,(-self.slit_size/2),(self.slit_size/2))
		self.percent_l,self.percent_err_l = integrate.quad(funx,(-self.seeing/2),(self.seeing/2))
		self.percent = self.percent_u * self.percent_l # can use error if you add it later...
		self.extension = self.seeing * self.slit_size
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

	def recalculate_dichroic(self,caller):
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
		
		rn = edl.rn_default

		if (self.bin_size > 0) and (self.bin_size < 5):
			self.read_noise = math.ceil(rn * spectral_resolution * spatial_resolution / (self.bin_size**2))
			if not self.sss:
				print('{} Pixel binning: ({}x{})'.format(edl.string_prefix,self.bin_size,self.bin_size))
				print('{} Extent: {} arcsec^2\n{} num pixels/resel: {} px\n{} spectral resolution: {} px\n{} spatial resolution: {} px'.format(edl.string_prefix,extent,edl.string_prefix,int(math.ceil(npix)),edl.string_prefix,spectral_resolution,edl.string_prefix,spatial_resolution))
		else:
			raise ValueError('{} Invalid pixel binning option: {}'.format(edl.string_prefix,self.bin_size))

	def recalculate_atmospheric_extinction(self):
		self.extinction = spectres(self.wavelength,dh.atmo_ext_x,dh.atmo_ext_y)