# this uses a caller singleton and enumerated caller-association groups to update only necessary data.
# the else statement catches errors and changes from wavelength or channel

import math
import numpy as np
import values as edl
import datahandler as dh
from spectres import spectres
from scipy import interpolate, integrate
from astropy import constants as const
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_sigma_to_fwhm

# handles various input forms
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
		bin_option=edl.bin_options_int[edl.bin_options_default_index],channel='both',sss=True):
		
		object_type = type_router(object_type)
		grating_opt = type_router(grating_opt)
		telescope_mode = type_router(telescope_mode)		
		self.plot_step = wavelength[2] - wavelength[1]
		# run init


	''' tier 0 '''

	def refresh(self,caller=cb_obj.name,tabs=tabs):
		# direct the changed values
		if (tabs.active in [_opt for _opt in range(3)]):
			if caller in edl.signal_components:
				self.recalculate_signal(caller)
			if caller in edl.noise_components:
				self.recalculate_noise(caller)
			if caller in edl.error_components:
				self.recalculate_error(caller)
			else:
				self.recalculate_signal(caller)
				self.recalculate_noise(caller)
				self.recalculate_error(caller)
		
			if (tabs.active == 0): # snr
				self.recalculate_snr(caller)
				plot_y1 = seflf.snr_blue
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
				raise ValueError("{} Invalid active tab: {}".format(string_prefix,tabs.active))

		return plot_y1,plot_y2


	def recalculate_snr(self,caller):
		if (self.channel == 'blue') or (self.channel == 'both'):
			self.snr_blue = np.divide(self.signal_blue,np.sqrt(self.signal_blue + self.noise_blue + np.square(self.readnoise)))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.snr_red = np.divide(self.signal_red,np.sqrt(self.signal_red + self.noise_red + np.square(self.readnoise)))


	def recalculate_error(self,caller):
		if (self.self.channel == 'blue') or (self.channel == 'both'):
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
		if caller in power_keys:
			self.recalculate_flux(caller)
			self.change_telescope_mode(caller)
		else:
			self.recalculate_flux(caller)
			self.change_telescope_mode(caller)

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
		if caller in edl.atmo_ext_keys:
			self.recalculate_atmospheric_extinction(caller)
		else:
			self.recalculate_grating(caller)
			self.recalculate_dichroic(caller)
			self.recalculate_ccd(caller)
			self.recalculate_atmospheric_extinction(caller)

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
		else:
			self.recalculate_dichroic(caller)
			self.recalculate_grating(caller)
			self.recalculate_ccd(caller)

		if (self.channel == 'blue') or (self.channel == 'both'):
			self.total_eff_noise_red = np.multiply(np.multiply(self.dichro_blue,self.grating_blue),(self.ccd_blue * np.square(self.mirror) * edl.coating_eff_blue))
		if (self.channel == 'red') or (self.channel == 'both'):
			self.total_eff_noise_blue = np.multiply(np.multiply(self.dichro_red,self.grating_red),(self.ccd_red * np.square(self.mirror) * edl.coating_eff_red))


	def recalculate_flux(self,caller):
		if caller in edl.grating_opt_keys:
			self.change_grating_opt(caller)
		if caller in edl.filter_keys:
			self.change_filter(caller)
		if caller in edl.mag_sys_keys():
			self.change_mag_sys_opt(caller)
		if caller in edl.atmo_ext_keys:
			self.change_moon_days(caller)
		else:
			self.change_grating_opt()
			self.change_filter()
			self.change_mag_sys_opt()
			self.change_moon_days()

		# heal identicalities
		self.lambda_A[0] = lambda_A[0] + self.plot_step
		self.lambda_A[-1] = lambda_A[-1] - self.plot_step

		ftrans = interpolate.interp1d(selected_filter[0],selected_filter[1], kind='cubic')
		trans = ftrans(self.lambda_A)

		_extinction = spectres(self.lambda_A,dh.atmo_ext_x,dh.atmo_ext_y)

		flux = self.flux_A * 1e10
		_lambda = self.lambda_A / 1e10

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
			raise ValueError("{} Invalid grating_opt: {}".format(string_prefix,self.grating_opt))
		if not sss:
			print("{} Grating changed to {}".format(string_prefix,_active))


	def change_mag_sys_opt(self,caller):
		if (self.mag_sys_opt == 'vega'):
			flux_vega = spectres(self.wavelength,dh.vega[0],dh.vega[1]) * 1e10 # fixed... I hope?
			self.mag_model = -2.5 * np.log10(np.divide(math.fsum(flux * _extinction * _lambda * trans),math.fsum(flux_vega * trans * _lambda * _extinction))) + 0.03
		elif (self.mag_sys_opt == 'ab'):
			self.mag_model = -48.6 - 2.5 * np.log10(math.fsum(flux * trans * _extinction *_lambda) / math.fsum(trans * _lambda * _extinction * (const.c.value/np.square(_lambda))))
		else:
			raise ValueError("{} Invalid magnitude system option (mag_sys_opt): {}".format(string_prefix,self.mag_sys_opt))
		if not sss:
			print("{} Magnitude system changed to {}".format(string_prefix,self.mag_sys_opt.upper()))


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
		if not sss:
			print("{} Object type changed to {}".format(string_prefix,self.object_type))


	def change_moon_days(self,caller):
		if self.moon_days in moon_days_keys:
			self.sky_background = dh.skyfiles[(int(np.where(np.asarray(edl.moon_days_keys)==self.moon_days)[0]))]
		else:
			raise ValueError('{} Invalid number of days since new moon: {}'.format(string_prefix,self.moon_days))
		self.recalculate_sky_flux(caller)
		if not sss:
			print("{} Days since new moon changed to {}".format(string_prefix,self.moon_days))


	def change_telescope_mode(self,caller):
		if self.telescope_mode in edl.telescope_mode_keys[:3]:
			self.area = edl.area[0]
		elif self.telescope_mode in edl.telescope_mode_keys[4:]:
			self.area = edl.area[1]
		else:
			raise ValueError('{} Invalid telescope mode: {}'.format(string_prefix,self.telescope_mode))
		if not sss:
			print("{} Telescope mode changed to {} (area: {} m^2)".format(string_prefix,self.telescope_mode,self.area))


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
		if not sss:
			_active = edl.filter_files[self.filter_index]
			print("{} Filter changed to {} ({}--{} nm)".format(string_prefix,_active,_active[0],_active[-1]))


	def recalculate_seeing(self,caller):
		_sigma = self.seeing / gaussian_sigma_to_fwhm
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		self.percent_u,self.percent_err_u = integrate.quad(funx,(-slit_size/2),(slit_size/2))
		self.percent_l,self.percent_err_l = integrate.quad(funx,(-seeing/2),(seeing/2))
		self.percent = self.percent_u * self.percent_l # can use error if you add it later...
		self.extension = self.seeing * self.slit_size
		if not sss:
			print("{} Seeing changed to {}".format(string_prefix,self.seeing))


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


	def recalculate_readnoise(self,caller): # probably not implemented in all necessary places yet...
		spectral_resolution = math.ceil((self.slit_size/(0.7/12))/2)*2 # px (ceil()/2)*2 to round up to next even integer
		spatial_resolution = math.ceil((self.seeing/(0.7/12))/2)*2 # probably have to find another method, because pscript is depriciating :(
		extent = seeing * slit_size
		npix = extent/(0.7/12)**2
		
		rn = edl.rn_default

		if (self.bin_size > 0) and (self.bin_size < 5):
			self.readnoise = math.ceil(rn * spectral_resolution * spatial_resolution / (self.bin_size**2))
			if not sss:
				print('{} Pixel binning: ({}x{})'.format(string_prefix,self.bin_size,self.bin_size))
				print('{} Extent: {} arcsec^2\n{} num pixels/resel: {} px\n{} spectral resolution: {} px\n{} spatial resolution: {} px'.format(string_prefix,extent,string_prefix,int(math.ceil(npix)),string_prefix,spectral_resolution,string_prefix,spatial_resolution))
		else:
			raise ValueError('{} Invalid pixel binning option ({})'.format(string_prefix,self.bin_size))


	def recalculate_atmospheric_extinction(self,caller):
		self.extinction = spectres(self.wavelength,dh.atmo_ext_x,dh.atmo_ext_y)


if __name__ == '__main__':
	main()