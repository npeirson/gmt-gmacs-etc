# this uses a caller singleton and enumerated caller-association groups to update only necessary data.
# the else statement catches errors and changes from wavelength or channel

import math
import numpy as np
import values as edl
import defaults as dfs
import datahandler as dh
from spectres import spectres
from scipy import interpolate, integrate
from astropy import constants as const
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_sigma_to_fwhm


# handles various input forms
def num_router(_obj):
	try:
		if isinstance(_obj,int):
			return _obj
		elif isinstance(_obj,float):
			return _obj
		else:
			raise ValueError("{} Invalid value: {}".format(edl.string_prefix,_obj))
	except:
		raise TypeError("{} Invalid type: {}".format(edl.string_prefix,type(_obj)))



def type_router(_obj,_keys):
	try:
		if isinstance(_obj,str):
			return int(np.where(np.asarray(_keys)==_obj.lower())[0])
		else:
			return num_router(_obj)
	except:
		raise TypeError("{} Invalid type: {}".format(edl.string_prefix,type(_obj)))



class calculator:
	"""
	Class for simulating observatations using GMACS, designed to serve as the back-end to a static bokeh site.

	Parameters
	----------


	"""

	def __init__(self,mode='snr',telescope_mode='first',wavelength=dfs.default_wavelength,exposure_time=3600,object_type='a5v',
		filter_index=3,mag_sys_opt='ab',magnitude=25,redshift=0,seeing=0.5,slit_size=0.5,moon_days=0,grating_opt=0,
		noise=False,bin_option=edl.bin_options_int[edl.bin_options_default_index],channel='both',sss=True,**kwargs):
		
		if isinstance(mag_sys_opt,str):
			if (mag_sys_opt.lower() == 'ab') or (mag_sys_opt.lower() == 'vega'):
				self.mag_sys_opt = mag_sys_opt
		elif isinstance(mag_sys_opt,int):
			if (mag_sys_opt == 0):
				self.mag_sys_opt = 'ab'
			elif (mag_sys_opt == 1):
				self.mag_sys_opt = 'vega'
			else:
				raise ValueError('{} Invalid mag_sys_opt value: {}'.format(edl.string_prefix,mag_sys_opt))
		else:
			raise TypeError('{} Invalid mag_sys_opt value type!'.format(edl.string_prefix))



		self.sss = sss
		self.mode = mode
		# required args, maybe add range limit too
		self.object_type = type_router(object_type,edl.object_type_keys)
		if (self.object_type >= len(edl.object_type_keys)):
			self.object_type = self.object_type - len(edl.object_type_keys)

		self.grating_opt = type_router(grating_opt,edl.grating_opt_keys)
		if (self.grating_opt >= 3):
			self.grating_opt = 1
		else:
			self.grating_opt = 0

		self.telescope_mode = type_router(telescope_mode,edl.telescope_mode_keys)
		if (self.telescope_mode >= 3):
			edl.telescope_mode = 1
		else:
			edl.telescope_mode = 0

		self.filter_index = type_router(filter_index,edl.filter_keys)
		self.magnitude = num_router(magnitude)
		self.seeing = num_router(seeing)
		self.slit_size = num_router(slit_size)
		self.redshift = num_router(redshift)
		self.moon_days = num_router(moon_days)
		if isinstance(wavelength,np.ndarray):
			self.wavelength = wavelength			
			if (self.wavelength[0] >= dfs.default_limits_wavelength[0]) and (self.wavelength[-1] <= dfs.default_limits_wavelength[1]):
				self.plot_step = self.wavelength[2] - self.wavelength[1]
			else:
				raise ValueError('{} Invalid wavelength extrema ({}--{}). Must be within {}--{}'.format(edl.string_prefix,wavelength[0],wavelength[-1],dfs.default_limits_wavelength[0],dfs.default_limits_wavelength[1]))
		else:
			raise TypeError("{} Invalid wavelength type: {} (must be array-like)".format(edl.string_prefix,type(wavelength)))
		#self.__dict__ = dict([arg for arg in kwargs if arg in edl.keys]) # pick up any optionals
		self.change('wavelength') # initializing plots is as easy as this


	''' tier 0 '''

	def change(self,caller): # caller=cb_obj.name,tabs=tabs
		# direct the changed values
		if (self.mode in edl.mode_keys.flatten()):
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
		
			if self.mode in edl.mode_keys[0]: # snr
				self.recalculate_snr(caller)
				plot_y1 = self.snr_blue
				plot_y2 = self.snr_red
			elif self.mode in edl.mode_keys[1]: # obs spec no noise
				if self.noise:
					plot_y1 = np.add(self.signal_blue,self.error_blue)
					plot_y2 = np.add(self.signal_red,self.error_blue)
				else:
					plot_y1 = self.signal_blue
					plot_y2 = self.signal_red
			elif self.mode in edl.mode_keys[2]: # obs sky background
				plot_y1 = self.error_blue
				plot_y2 = self.error_red
			else:
				pass # pass to other options

		else:
			if self.mode in edl.mode_keys[3]: # dichroic throughput
				plot_y1 = self.dichro_blue
				plot_y2 = self.dichro_red
			elif self.mode in edl.mode_keys[4]: # grating throughput
				plot_y1 = self.grating_blue
				plot_y2 = self.grating_red
			elif self.mode in edl.mode_keys[5]: # CCD QE
				plot_y1 = self.ccd_blue
				plot_y2 = self.ccd_red
			elif self.mode in edl.mode_keys[6]: # atmospheric extinction
				plot_y1 = self.extinction
				plot_y2 = self.extinction
			else:
				raise ValueError("{} Invalid active tab: {}".format(edl.string_prefix,self.mode))

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
		else:
			self.recalculate_signal(caller)
			self.recalculate_noise(caller)
			self.recalculate_readnoise(caller)

		if (self.self.channel == 'blue') or (self.channel == 'both'):
			sigma_blue = np.sqrt(self.signal_blue + self.noise_blue + np.square(self.readnoise))
			self.error_blue = np.random.normal(loc=0, scale=sigma_blue,size=len(self.snr_blue))
		if (self.channel == 'red') or (self.channel == 'both'):
			sigma_red = np.sqrt(self.signal_red + self.noise_red + np.square(self.readnoise))
			self.error_red = np.random.normal(loc=0, scale=sigma_red,size=len(self.snr_red))


	''' tier 1 '''

	def recalculate_signal(self,caller):
		if caller in edl.counts_keys:
			self.recalculate_counts(caller)
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
		elif caller in edl.total_eff_noise_keys:
			self.recalculalte_efficiency_noise(caller)
		else:
			self.recalculate_counts_noise(caller)
			self.recalculalte_efficiency_noise(caller)

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.noise_blue = np.multiply(self.counts_noise,self.total_eff_noise_blue)
		if (channel == 'red') or (self.channel == 'both'):
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
			self.recalculate_extension(caller)
		else:
			self.recalculate_sky_flux(caller)
			self.recalculate_extension(caller)
			self.change_telescope_mode(caller)
			self.change_plot_step(caller)

		self.counts_noise = np.multiply(np.multiply(self.sky_flux,self.extension),(self.area*self.exposure_time*self.plot_step))


	def recalculate_percent(self,caller):
		self.recalculate_seeing(caller)
		# this is a forwarding function


	def recalculate_efficiency(self,caller):
		if caller in edl.grating_keys:
			self.recalculate_grating(caller)
		if caller in edl.dichroic_keys:
			self.recalculate_dichroic(caller)
		if caller in edl.ccd_keys:
			self.recalculate_ccd(caller)
		if caller in edl.mirror_keys:
			self.recalculate_mirror(caller)
		else:
			self.recalculate_grating(caller)
			self.recalculate_dichroic(caller)
			self.recalculate_ccd(caller)
			self.recalculate_atmospheric_extinction(caller)
			self.recalculate_mirror(caller)

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.total_eff_blue = np.multiply(np.multiply(self.dichro_blue,self.grating_blue),np.multiply((self.ccd_blue * (edl.coating_eff_blue * self.extinction)),np.square(mirror)))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.total_eff_red = np.multiply(np.multiply(self.dichro_red,self.grating_red),np.multiply((self.ccd_red * (edl.coating_eff_red * self.extinction)),np.square(mirror)))


	def recalculalte_efficiency_noise(self,caller):
		if caller in edl.dichroic_keys:
			self.recalculate_dichroic(caller)
		if caller in edl.grating_keys:
			self.recalculate_grating(caller)
		if caller in edl.ccd_keys:
			self.recalculate_ccd(caller)
		if caller in edl.mirror_keys:
			self.recalculate_mirror(caller)
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
		if caller in edl.object_type_keys:
			self.change_object_type(caller)
		else:
			self.change_object_type(caller)
			self.change_grating_opt(caller)
			self.change_filter(caller)
			self.change_moon_days(caller)
			self.change_plot_step(caller)

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
		if (self.mag_sys_opt == 'vega'):
			flux_vega = spectres(self.wavelength,dh.vega[0],dh.vega[1]) * 1e10 # fixed... I hope?
			self.mag_model = -2.5 * np.log10(np.divide(math.fsum(self.flux * self._extinction * self._lambda * self.trans),math.fsum(self.flux_vega * self.trans * self._lambda * self._extinction))) + 0.03
		elif (self.mag_sys_opt == 'ab'):
			self.mag_model = -48.6 - 2.5 * np.log10(math.fsum(self.flux * self.trans * self._extinction * self._lambda) / math.fsum(self.trans * self._lambda * self._extinction * (const.c.value/np.square(self._lambda))))
		else:
			raise ValueError("{} Invalid magnitude system option (mag_sys_opt): {}".format(edl.string_prefix,self.mag_sys_opt))
		if not self.sss:
			print("{} Magnitude system changed to {}".format(edl.string_prefix,self.mag_sys_opt.upper()))


	def change_object_type(self,caller):
		if self.object_type in edl.stellar_keys:
			index_of = [i for i,name in enumerate(stellar_keys) if self.object_type in name][0]
			self.object_type = dh.starfiles[index_of]
		elif self.object_type in edl.galactic_keys:
			index_of = [i for i,name in enumerate(galactic_keys) if self.object_type in name][0]
			self.object_type = dh.galaxyfiles[index_of]
		else:
			raise ValueError("{} Invalid object type: {}".format(edl.string_prefix,self.object_type))

		self.object_x = self.object_type[0] * (1+redshift)
		self.object_y = self.object_type[1]
		self.flux_A = spectres(lambda_A,object_x,object_y)
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
		lambda_min = filter_min
		filter_max = max(self.selected_filter[0])
		lambda_max = filter_max

		self.change_plot_step(caller)
		self.lambda_A = np.arange(lambda_min,lambda_max,self.plot_step)
		if not self.sss:
			_active = edl.filter_files[self.filter_index]
			print("{} Filter changed to {} ({}--{} nm)".format(edl.string_prefix,_active,self.selected_filter[0],self.selected_filter[-1]))


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


	def recalculate_grating(self,caller):
		if (self.channel is 'blue') or (self.channel is 'both'):
			self.grating_blue = spectres(self.wavelength,(dh.grating1[0]*10),dh.grating1[1])
		if (self.channel is 'red') or (self.channel is 'both'):
			self.grating_red = spectres(self.wavelength,(dh.grating2[0]*10),dh.grating2[1])


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
		extent = seeing * slit_size
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
		self.extinction = spectres(self.wavelength,dh.atmo_ext_x,dh.atmo_ext_y)


	def change_plot_step(self,caller):
		if (self.wavelength[0] >= dfs.default_limits_wavelength[0]) and (self.wavelength[-1] <= dfs.default_limits_wavelength[1]):
			self.plot_step = self.wavelength[2] - self.wavelength[1]
		else: # technically shouldn't happen...
			raise ValueError('{} Invalid wavelength extrema ({}--{}). Must be within {}--{}'.format(edl.string_prefix,wavelength[0],wavelength[-1],dfs.default_limits_wavelength[0],dfs.default_limits_wavelength[1]))