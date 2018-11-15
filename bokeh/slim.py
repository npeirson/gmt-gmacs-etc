#!/usr/bin/python3

''' GMACS Exposure Time Calculator SLIM '''

import math
import time
import numpy as np

from spectres import spectres
from scipy import interpolate, integrate
from astropy import constants # `const` not permitted in JS
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_sigma_to_fwhm
from bokeh.models import ColumnDataSource

import datahandler as dh
import strings as stc
import tools as etools
import keys as etkeys
import paths as etpaths
import defaults as dfs


class session:

	def __init__(self,*args):
		input_dict = dict()
		print(range(len(args)))
		for i in range(len(args)):
			input_dict.update(args[i])
		self.__dict__ = input_dict

		# initialize values
		self.sss = True # False for verbose mode
		self.wavelength_array = np.arange(self.wavelength.value[0],self.wavelength.value[1],dfs.plot_step)
		self.init()

	def init(self):
		self.time_last_command = time.time()
		return self.snr('init') # initialize

		
	def update(self,caller):
		if ((time.time() - self.time_last_command) >= .01):
			self.time_last_command = time.time()
			print("{} Call from: {}".format(dfs.string_prefix,caller))

			if (self.tabs.active == 0):
				return self.snr(caller)
			elif (self.tabs.active == 1):
				return self.os(caller)
			elif (self.tabs.active == 2):
				return self.sky(caller)
			elif (self.tabs.active == 3):
				return self.dichro(caller)
			elif (self.tabs.active == 4):
				return self.gt(caller)
			elif (self.tabs.active == 5):
				return self.ccd(caller)
			elif (self.tabs.active == 6):
				return self.atmo_ext(caller)
			else:
				raise ValueError("{} Invalid update request!".format(dfs.string_prefix))
		else:
			if not self.sss:
				print('{} Caught Buffer overflow!'.format(dfs.string_prefix))


	''' mode functions '''
	def snr(self,caller):	
		x = self.wavelength_array
		yb,yr = self.calc_snr(caller)
		return x,yb,yr

	def os(self,caller):
		x = self.wavelength_array
		signal_blue,signal_red = self.recalculate_signal(caller)
		if self.withnoise.active:
			error_blue,error_red = self.calc_error(caller)
			yb,yr = np.add(signal_blue,error_blue),np.add(signal_red,error_red)
		else:
			yb,yr = signal_blue,signal_red
		return x,yb,yr

	def sky(self,caller):
		x = self.wavelength_array
		yb,yr = self.calc_error(caller)
		return x,yb,yr

	def dichro(self,caller):
		x = self.wavelength_array
		yb,yr = self.recalculate_dichroic(caller)
		return x,yb,yr

	def gt(self,caller):
		x = self.wavelength_array
		yb,yr = self.recalculate_grating(caller)
		return x,yb,yr

	def ccd(self,caller):
		x = self.wavelength_array
		yb,yr = self.recalculate_ccd(caller)
		return x,yb,yr

	def atmo_ext(self,caller):
		x = self.wavelength_array
		yb,yr = self.recalculate_atmospheric_extinction(caller)
		return x,yb,yr


	def calc_snr(self,caller):
		read_noise = self.recalculate_readnoise(caller)
		signal_blue,signal_red = self.recalculate_signal(caller)
		noise_blue,noise_red = self.recalculate_noise(caller)
		snr_blue = np.divide(signal_blue,np.sqrt(signal_blue + noise_blue + np.square(read_noise)))
		snr_red = np.divide(signal_red,np.sqrt(signal_red + noise_red + np.square(read_noise)))
		return snr_blue,snr_red

	def calc_error(self,caller):
		snr_blue,snr_red = self.calc_snr(caller)
		read_noise = self.recalculate_readnoise(caller)
		signal_blue,signal_red = self.recalculate_signal(caller)
		noise_blue,noise_red = self.recalculate_noise(caller)
		sigma_blue = np.sqrt(signal_blue + noise_blue + np.square(read_noise))
		error_blue = np.random.normal(loc=0, scale=sigma_blue,size=len(snr_blue))
		sigma_red = np.sqrt(signal_red + noise_red + np.square(read_noise))
		error_red = np.random.normal(loc=0, scale=sigma_red,size=len(snr_red))
		return error_blue,error_red


	def recalculate_signal(self,caller):
		counts = self.recalculate_counts(caller)
		percent,extension = self.recalculate_seeing(caller)
		total_eff_blue,total_eff_red = self.recalculate_efficiency(caller)
		signal_blue = np.multiply((counts * percent), total_eff_blue)
		signal_red = np.multiply((counts * percent), total_eff_red)
		return signal_blue,signal_red

	def recalculate_noise(self,caller):
		counts_noise = self.recalculate_counts_noise(caller)
		total_eff_noise_blue,total_eff_noise_red = self.recalculate_efficiency_noise(caller)
		noise_blue = np.multiply(counts_noise,total_eff_noise_blue)
		noise_red = np.multiply(counts_noise,total_eff_noise_red)
		return noise_blue,noise_red

	def recalculate_counts(self,caller):
		self.recalculate_wavelength(caller)
		flux,flux_y = self.recalculate_flux(caller)
		area = self.change_num_mirrors(caller)
		power = flux_y * area * self.time.value * self.plot_step
		counts = np.divide(np.divide(power,np.divide((constants.h.value * constants.c.value),self.wavelength_array)),1e10)
		return counts

	def recalculate_counts_noise(self,caller):
		self.recalculate_wavelength(caller)
		sky_flux = self.recalculate_sky_flux(caller)
		percent,extension = self.recalculate_seeing(caller)
		area = self.change_num_mirrors(caller)
		counts_noise = np.multiply(np.multiply(sky_flux,extension),(area*self.time.value*self.plot_step))
		return counts_noise

	def recalculate_efficiency(self,caller):
		grating_blue,grating_red = self.recalculate_grating(caller)
		dichro_blue,dichro_red = self.recalculate_dichroic(caller)
		ccd_blue,ccd_red = self.recalculate_ccd(caller)
		extinction = self.recalculate_atmospheric_extinction(caller)
		mirror_loss = self.recalculate_mirror_loss(caller)
		total_eff_blue = np.multiply(np.multiply(dichro_blue,grating_blue),np.multiply((ccd_blue * (dfs.coating_eff_blue * extinction)),np.square(mirror_loss)))
		total_eff_red = np.multiply(np.multiply(dichro_red,grating_red),np.multiply((ccd_red * (dfs.coating_eff_red * extinction)),np.square(mirror_loss)))
		return total_eff_blue,total_eff_red

	def recalculate_efficiency_noise(self,caller):
		dichro_blue,dichro_red = self.recalculate_dichroic(caller)
		grating_blue,grating_red = self.recalculate_grating(caller)
		ccd_blue,ccd_red = self.recalculate_ccd(caller)
		mirror_loss = self.recalculate_mirror_loss(caller)
		total_eff_noise_red = np.multiply(np.multiply(dichro_blue,grating_blue),(ccd_blue * np.square(mirror_loss) * dfs.coating_eff_blue))
		total_eff_noise_blue = np.multiply(np.multiply(dichro_red,grating_red),(ccd_red * np.square(mirror_loss) * dfs.coating_eff_red))
		return total_eff_noise_blue,total_eff_noise_red

	def recalculate_flux(self,caller):
		self.recalculate_wavelength(caller)
		grating_blue,grating_red = self.recalculate_grating(caller)
		moon_days = self.change_moon_days(caller)
		selected_filter_x,selected_filter_y,lambda_A = self.change_filter(caller)
		object_x,object_y,flux_A = self.change_object_type(caller)

		# heal identicalities
		lambda_A[0] = lambda_A[0] + self.plot_step
		lambda_A[-1] = lambda_A[-1] - self.plot_step

		# recalculate some losses
		ftrans = interpolate.interp1d(selected_filter_x,selected_filter_y, kind='cubic')
		trans = ftrans(lambda_A)
		_extinction = spectres(lambda_A,dh.atmo_ext_x,dh.atmo_ext_y)
		_lambda = lambda_A / 1e10
		flux = flux_A * 1e10
		# valaidate flux
		num_zeros = 0
		for lux in flux:
			if (lux is None) or (lux is 0):
				num_zeros += 1
				lux = 0
		if (num_zeros >= (flux.shape[0]/5)):
			if (num_zeros == flux.shape[0]):
				print('{} No flux in this bandpass!'.format(dfs.string_prefix))
				output_flux = np.zeros(self.wavelength_array.shape[0])
			else:
				percent_zeros = (num_zeros / flux.shape[0]) * 100
				print('{} {}% of this bandpass has zero flux'.format(dfs.string_prefix,percent_zeros))
		# filter bandpass and associated flux
		if (self.mag_sys.active == 0):
			flux_vega = spectres(self.wavelength_array,dh.vega_file[0],dh.vega_file[1]) * 1e10 # fixed... I hope?
			mag_model = -2.5 * np.log10(np.divide(math.fsum(flux * _extinction * _lambda * trans),math.fsum(flux_vega * trans * _lambda * _extinction))) + 0.03
		elif (self.mag_sys.active == 1):
			mag_model = -48.6 - 2.5 * np.log10(math.fsum(flux * trans * _extinction * _lambda) / math.fsum(trans * _lambda * _extinction * (constants.c.value/np.square(_lambda))))
		else:
			raise ValueError("{} Invalid magnitude system option (mag_sys_opt): {}".format(dfs.string_prefix,self.mag_sys.active))	
		del_mag = self.mag.value - mag_model
		output_flux = np.multiply(object_y,10 ** np.negative(del_mag/2.5))
		old_res = object_x[2] - object_x[1]
		if (old_res < self.plot_step):
			flux_y = spectres(self.wavelength_array,object_x,(output_flux*1e-03)) # ergs s-1 cm-2 A-1 to J s-1 m-2 A-1
		else:
			flux_y = spectres(self.wavelength_array,object_x,(output_flux*1e-03))
		return flux,flux_y
		
	def recalculate_wavelength(self,caller):
		if (caller == 'init') or (caller == 'wavelength'):
			self.wavelength_array = np.arange(self.wavelength.value[0],self.wavelength.value[1],dfs.plot_step) # change plot step later if dynamic algorithm desired
			self.plot_step = self.wavelength_array[2] - self.wavelength_array[1]
		else:
			print('arg')

	def change_grating_opt(self,caller):
		delta_lambda = dfs.dld[self.grating.active] * self.slit.value / 0.7
		if not self.sss:
			print("{} Grating: \'{} resolution\'".format(dfs.string_prefix,'high' if (self.grating.active==1) else 'low'))
		return delta_lambda

	def change_object_type(self,caller):
		if (self.star_type.value == None):
			curr_star,curr_gal = 4,0
		else:
			curr_star = int(np.where(np.asarray(stc.star_types)==self.star_type.value)[0])
			curr_gal = int(np.where(np.asarray(stc.galaxy_types)==self.galaxy_type.value)[0])
		selected_filter_x,selected_filter_y,lambda_A = self.change_filter(caller)
		if (self.object_type.active == 0): # stellar classifications
			object_data = dh.starfiles[curr_star]
		elif (self.object_type.active == 1):
			object_data = dh.starfiles[curr_gal]
		else:
			raise ValueError("{} Invalid object class: ({})".format(dfs.string_prefix,stc.object_types[self.object_type.active]))
		object_x = object_data[0] * (1 + self.redshift.value)
		object_y = object_data[1]
		flux_A = spectres(lambda_A,object_x,object_y)
		if not self.sss:
			print("{} Object type: {}".format(dfs.string_prefix,stc.object_types[self.object_type.active]))
		return object_x,object_y,flux_A
			
	def change_moon_days(self,caller):
		if not self.sss:
			print("{} Days since/until new moon: {}".format(dfs.string_prefix,stc.moon_opts[self.moon.active]))
		return self.recalculate_sky_flux(caller)

	def change_num_mirrors(self,caller):
		if (self.telescope.active == 0):
			return dfs.area[0]
		else:
			return dfs.area[1]

	def change_filter(self,caller):
		if (self.filter.value == None):
			curr_fil = 7
		else:
			curr_fil = int(np.where(np.asarray(etpaths.filter_files)==self.filter.value)[0])
		selected_filter = dh.filterfiles[curr_fil]
		filter_min = min(selected_filter[0])
		filter_max = max(selected_filter[0])
		lambda_min,lambda_max = filter_min,filter_max
		self.recalculate_wavelength(caller)
		lambda_A = np.arange(lambda_min,lambda_max,self.plot_step)
		if not self.sss:
			_active = stc.filter_opts[curr_fil[1]]
			print("{} Filter: {} ({}--{} nm)".format(dfs.string_prefix,stc.filter_opts[curr_fil[0]],selected_filter_x[0],selected_filter_x[-1]))
		return selected_filter[0],selected_filter[1],lambda_A


	def recalculate_seeing(self,caller):
		_sigma = self.seeing.value / gaussian_sigma_to_fwhm
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		percent_u,percent_err_u = integrate.quad(funx,(-self.slit.value/2),(self.slit.value/2))
		percent_l,percent_err_l = integrate.quad(funx,(-self.seeing.value/2),(self.seeing.value/2))
		percent = percent_u * percent_l # can use error if you add it later...
		extension = self.seeing.value * self.slit.value
		return percent,extension

	def recalculate_sky_flux(self,caller):		
		self.recalculate_wavelength(caller)
		delta_lambda = self.change_grating_opt(caller)
		sky_x = dh.skyfiles[self.moon.active][0] * 1e4
		sky_y = dh.skyfiles[self.moon.active][1] / 1e4
		old_res = sky_x[2] - sky_x[1]
		_sigma = delta_lambda / gaussian_sigma_to_fwhm
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		_x = np.arange((-5*_sigma),(5*_sigma),old_res)
		degrade = funx(_x)/np.trapz(funx(_x))
		sky_y = convolve_fft(sky_y,degrade)
		return spectres(self.wavelength_array,sky_x,sky_y)

	def recalculate_dichroic(self,caller):
		self.recalculate_wavelength(caller)
		fblue_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y1, kind='cubic')
		fred_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y2, kind='cubic')
		return fblue_dichro(self.wavelength_array),fred_dichro(self.wavelength_array)

	def recalculate_grating(self,caller):
		self.recalculate_wavelength(caller)
		grating_blue = spectres(self.wavelength_array,(dh.grating1[0]*10),dh.grating1[1])
		grating_red = spectres(self.wavelength_array,(dh.grating2[0]*10),dh.grating2[1])
		return grating_blue,grating_red

	def recalculate_ccd(self,caller):
		self.recalculate_wavelength(caller)
		fblue_ccd = interpolate.interp1d((dh.ccd1[0]*10),dh.ccd1[1], kind='cubic')
		fred_ccd = interpolate.interp1d((dh.ccd2[0]*10),dh.ccd2[1], kind='cubic')
		return fblue_ccd(self.wavelength_array),fred_ccd(self.wavelength_array)

	def recalculate_mirror_loss(self,caller):
		self.recalculate_wavelength(caller)
		fmirror = interpolate.interp1d(dh.mirror_file_x,dh.mirror_file_y, kind='cubic')
		return fmirror(self.wavelength_array)

	def recalculate_readnoise(self,caller): # probably not implemented in all necessary places yet...
		spectral_resolution = math.ceil((self.slit.value/(0.7/12))/2)*2 # px (ceil()/2)*2 to round up to next even integer
		spatial_resolution = math.ceil((self.seeing.value/(0.7/12))/2)*2 # probably have to find another method, because pscript is depriciating :(
		extent = self.seeing.value * self.slit.value
		npix = extent/(0.7/12)**2
		rn = dfs.rn_default
		read_noise = math.ceil(rn * spectral_resolution * spatial_resolution / ((self.binning.active + 1)**2)) # +1 because zero indexing
		if not self.sss:
			print('{0} Pixel binning: ({1}x{1})'.format(dfs.string_prefix,(self.binning.active + 1)))
			print('{0} Extent: {1} arcsec^2\n{0} num pixels/resel: {2} px\n{0} spectral resolution: {3} px\n{0} spatial resolution: {4} px'.format(dfs.string_prefix,extent,int(math.ceil(npix)),spectral_resolution,spatial_resolution))
		return read_noise

	def recalculate_atmospheric_extinction(self,caller):
		extinction = spectres(self.wavelength_array,dh.atmo_ext_x,dh.atmo_ext_y)
		return extinction