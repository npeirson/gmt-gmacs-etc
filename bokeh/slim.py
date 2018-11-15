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
		for i in range(len(args)):
			input_dict.update(args[i])
		self.__dict__ = input_dict
		self.sss = True # False for verbose mode
		self.wavelength_array = np.arange(self.wavelength.value[0],self.wavelength.value[1],dfs.plot_step)
		self.plot_step = self.wavelength_array[2] - self.wavelength_array[1]
		self.recalculate_wavelength()
		# initialize columndatasources
		self.init()

	def init(self):
		return self.snr('init') # initialize

		
	def update(self,caller):
		print("{} Call from: {}".format(dfs.string_prefix,caller))
		if (caller == self.wavelength.name):
			self.recalculate_wavelength()
		if (self.tabs.active == 0):
			return self.snr(caller)
		elif (self.tabs.active == 1):
			return self.os()
		elif (self.tabs.active == 2):
			return self.sky()
		elif (self.tabs.active == 3):
			return self.dichro()
		elif (self.tabs.active == 4):
			return self.gt()
		elif (self.tabs.active == 5):
			return self.ccd()
		elif (self.tabs.active == 6):
			return self.atmo_ext()
		else:
			raise ValueError("{} Invalid update request!".format(dfs.string_prefix))


	'''
		mode functions
			- first collect needed values
			- then do maths
			this avoids unnecessary work

			and is also a complete mess

	'''

	def snr(self,caller):
		read_noise = self.recalculate_readnoise(caller)
		area = self.change_num_mirrors()
		percent,extension = self.recalculate_seeing()
		selected_filter_x,selected_filter_y,lambda_A = self.change_filter(caller)
		moon_days = self.change_moon_days()
		percent,extension = self.recalculate_seeing()
		dichro_blue,dichro_red = self.recalculate_dichroic(caller)
		grating_blue,grating_red = self.recalculate_grating(caller)
		ccd_blue,ccd_red = self.recalculate_ccd(caller)
		mirror_loss = self.recalculate_mirror_loss(caller)
		extinction = self.recalculate_atmospheric_extinction(caller)
		counts_noise = self.recalculate_counts_noise(sky_y=moon_days,extension=extension,area=area)

		total_eff_noise_blue,total_eff_noise_red = self.recalculate_efficiency_noise(dichro_blue=dichro_blue,dichro_red=dichro_red,grating_blue=grating_blue,grating_red=grating_red,ccd_blue=ccd_blue,ccd_red=ccd_red,mirror_loss=mirror_loss)
		noise_blue,noise_red = self.recalculate_noise(counts_noise=counts_noise,total_eff_noise_blue=total_eff_noise_blue,total_eff_noise_red=total_eff_noise_red)
		object_x,object_y,flux_A = self.change_object_type(selected_filter_x=selected_filter_x,selected_filter_y=selected_filter_y,lambda_A=lambda_A)
		flux,flux_y = self.recalculate_flux(grating_blue=grating_blue,grating_red=grating_red,moon_days=moon_days,selected_filter_x=selected_filter_x,selected_filter_y=selected_filter_y,lambda_A=lambda_A,object_x=object_x,object_y=object_y,flux_A=flux_A)
		counts = self.recalculate_counts(flux_y=flux_y,area=area)
		total_eff_blue,total_eff_red = self.recalculate_efficiency(grating_blue=grating_blue,grating_red=grating_red,dichro_blue=dichro_blue,dichro_red=dichro_red,ccd_blue=ccd_blue,ccd_red=ccd_red,extinction=extinction,mirror_loss=mirror_loss)
		signal_blue,signal_red = self.recalculate_signal(counts=counts,percent=percent,total_eff_blue=total_eff_blue,total_eff_red=total_eff_red)

		x = self.wavelength_array
		yb,yr = self.calc_snr(read_noise=read_noise,signal_blue=signal_blue,signal_red=signal_red,noise_blue=noise_blue,noise_red=noise_red)
		return x,yb,yr

	def os(self):
		object_x,object_y,flux_A = self.change_object_type(selected_filter_x=selected_filter_x,selected_filter_y=selected_filter_y,lambda_A=lambda_A)
		read_noise = self.recalculate_readnoise()
		area = self.change_num_mirrors()
		percent,extension = self.recalculate_seeing()
		dichro_blue,dichro_red = self.recalculate_dichroic()
		grating_blue,grating_red = self.recalculate_grating()
		ccd_blue,ccd_red = self.recalculate_ccd()
		mirror_loss = self.recalculate_mirror_loss()
		extinction = self.recalculate_atmospheric_extinction()
		selected_filter_x,selected_filter_y,lambda_A = self.change_filter()
		sky_y = self.change_moon_days()
		counts_noise = self.recalculate_counts_noise(sky_y=sky_y,extension=extension,area=area)
		noise_blue,noise_red = self.recalculate_noise(counts_noise=counts_noise,total_eff_noise_blue=total_eff_noise_blue,total_eff_noise_red=total_eff_noise_red)
		object_x,object_y,flux_A = self.change_object_type(selected_filter_x=selected_filter_x,selected_filter_y=selected_filter_y,lambda_A=lambda_A)
		flux,flux_y = self.recalculate_flux(grating_blue=grating_blue,grating_red=grating_red,moon_days=sky_y,selected_filter_x=selected_filter_x,selected_filter_y=selected_filter_y,lambda_A=lambda_A,object_x=object_x,object_y=object_y,flux_A=flux_A)
		counts = self.recalculate_counts(flux_y=flux_y,area=area)
		total_eff_blue,total_eff_red = self.recalculate_efficiency(grating_blue=grating_blue,grating_red=grating_red,dichro_blue=dichro_blue,dichro_red=dichro_red,ccd_blue=ccd_blue,ccd_red=ccd_red,extinction=extinction,mirror_loss=mirror_loss)
		signal_blue,signal_red = self.recalculate_signal(counts=counts,percent=percent,total_eff_blue=total_eff_blue,total_eff_red=total_eff_red)
		x = self.wavelength_array
		if self.withnoise.active:
			snr_blue,snr_red = self.calc_snr(read_noise=read_noise,signal_blue=signal_blue,signal_red=signal_red,noise_blue=noise_blue,noise_red=noise_red)
			error_blue,error_red = self.calc_error(snr_blue=snr_blue,snr_red=snr_red,read_noise=read_noise)
			yb,yr = np.add(signal_blue,error_blue),np.add(signal_red,error_red)
		else:
			yb,yr = signal_blue,signal_red
		return x,yb,yr

	def sky(self,caller):
		read_noise = self.recalculate_readnoise(caller)
		x,snr_blue,snr_red = self.snr()
		yb,yr = self.calc_error(snr_blue=snr_blue,snr_red=snr_red,read_noise=read_noise)
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


	def calc_snr(self,read_noise,signal_blue,signal_red,noise_blue,noise_red):
		snr_blue = np.divide(signal_blue,np.sqrt(signal_blue + noise_blue + np.square(read_noise)))
		snr_red = np.divide(signal_red,np.sqrt(signal_red + noise_red + np.square(read_noise)))
		return snr_blue,snr_red

	def calc_error(self,snr_blue,snr_red,read_noise):
		sigma_blue = np.sqrt(signal_blue + noise_blue + np.square(read_noise))
		error_blue = np.random.normal(loc=0, scale=sigma_blue,size=len(snr_blue))
		sigma_red = np.sqrt(signal_red + noise_red + np.square(read_noise))
		error_red = np.random.normal(loc=0, scale=sigma_red,size=len(snr_red))
		return error_blue,error_red


	def recalculate_signal(self,counts,percent,total_eff_blue,total_eff_red):
		signal_blue = np.multiply((counts * percent), total_eff_blue)
		signal_red = np.multiply((counts * percent), total_eff_red)
		return signal_blue,signal_red

	def recalculate_noise(self,counts_noise,total_eff_noise_blue,total_eff_noise_red):
		noise_blue = np.multiply(counts_noise,total_eff_noise_blue)
		noise_red = np.multiply(counts_noise,total_eff_noise_red)
		return noise_blue,noise_red

	def recalculate_counts(self,flux_y,area):
		power = flux_y * area * self.time.value * self.plot_step
		counts = np.divide(np.divide(power,np.divide((constants.h.value * constants.c.value),self.wavelength_array)),1e10)
		return counts

	def recalculate_counts_noise(self,sky_y,extension,area):
		counts_noise = np.multiply(np.multiply(sky_y,extension),(area*self.time.value*self.plot_step))
		return counts_noise

	def recalculate_efficiency(self,grating_blue,grating_red,dichro_blue,dichro_red,ccd_blue,ccd_red,extinction,mirror_loss):
		total_eff_blue = np.multiply(np.multiply(dichro_blue,grating_blue),np.multiply((ccd_blue * (dfs.coating_eff_blue * extinction)),np.square(mirror_loss)))
		total_eff_red = np.multiply(np.multiply(dichro_red,grating_red),np.multiply((ccd_red * (dfs.coating_eff_red * extinction)),np.square(mirror_loss)))
		return total_eff_blue,total_eff_red

	def recalculate_efficiency_noise(self,dichro_blue,dichro_red,grating_blue,grating_red,ccd_blue,ccd_red,mirror_loss):
		total_eff_noise_red = np.multiply(np.multiply(dichro_blue,grating_blue),(ccd_blue * np.square(mirror_loss) * dfs.coating_eff_blue))
		total_eff_noise_blue = np.multiply(np.multiply(dichro_red,grating_red),(ccd_red * np.square(mirror_loss) * dfs.coating_eff_red))
		return total_eff_noise_blue,total_eff_noise_red

	def recalculate_flux(self,grating_blue,grating_red,moon_days,selected_filter_x,selected_filter_y,lambda_A,object_x,object_y,flux_A):
		# heal identicalities, do this to main lams too
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
		
	def recalculate_wavelength(self):
		if (self.wavelength_array[0] != self.wavelength.value[0]) or (self.wavelength_array[-1] != self.wavelength.value[1]):
			self.wavelength_array = np.arange(self.wavelength.value[0],self.wavelength.value[1],dfs.plot_step) # change plot step later if dynamic algorithm desired
			self.plot_step = self.wavelength_array[2] - self.wavelength_array[1]

	def change_object_type(self,selected_filter_x,selected_filter_y,lambda_A):
		if (self.star_type.value == None):
			curr_star,curr_gal = 4,0
		else:
			curr_star = np.where(np.asarray(stc.star_types)==self.star_type.value)[0][0]
			curr_gal = np.where(np.asarray(stc.galaxy_types)==self.galaxy_type.value)[0][0]
		if (self.object_type.active == 0): # stellar classifications
			object_data = dh.starfiles[curr_star]
		elif (self.object_type.active == 1):
			object_data = dh.galaxyfiles[curr_gal]
		else:
			raise ValueError("{} Invalid object class: ({})".format(dfs.string_prefix,stc.object_types[self.object_type.active]))
		object_x = object_data[0] * (1 + self.redshift.value)
		object_y = object_data[1]
		flux_A = spectres(lambda_A,object_x,object_y)
		if not self.sss:
			print("{} Object type: {}".format(dfs.string_prefix,stc.object_types[self.object_type.active]))
		return object_x,object_y,flux_A


	def change_grating_opt(self,caller):
		if caller in etkeys.func_grating:
			self.stored_delta_lambda = dfs.dld[self.grating.active] * self.slit.value / 0.7
			if not self.sss:
				print("{} Grating: \'{} resolution\'".format(dfs.string_prefix,'high' if (self.grating.active==1) else 'low'))
		return self.stored_delta_lambda

			
	def change_moon_days(self):
		delta_lambda = self.change_grating_opt('init')
		if not self.sss:
			print("{} Days since/until new moon: {}".format(dfs.string_prefix,stc.moon_opts[self.moon.active]))
		return self.recalculate_sky_flux(delta_lambda)

	def change_num_mirrors(self):
		if (self.telescope.active == 0):
			return dfs.area[0]
		else:
			return dfs.area[1]

	def change_filter(self,caller):
		if (caller == 'init'):
			curr_fil = 7
		else:
			curr_fil = np.where(np.asarray(etpaths.filter_files)==self.filter.value)[0]
		if caller in etkeys.func_filter:
			selected_filter = dh.filterfiles[curr_fil]
			filter_min = min(selected_filter[0])
			filter_max = max(selected_filter[0])
			lambda_min,lambda_max = filter_min,filter_max
			lambda_A = np.arange(lambda_min,lambda_max,self.plot_step)
			self.store_filer_0,self.store_filter_1,self.store_filter_2 = selected_filter[0],selected_filter[1],lambda_A
			if not self.sss:
				_active = stc.filter_opts[curr_fil[1]]
				print("{} Filter: {} ({}--{} nm)".format(dfs.string_prefix,stc.filter_opts[curr_fil[0]],selected_filter_x[0],selected_filter_x[-1]))				
		return self.store_filer_0,self.store_filter_1,self.store_filter_2


	def recalculate_seeing(self):
		_sigma = self.seeing.value / gaussian_sigma_to_fwhm
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		percent_u,percent_err_u = integrate.quad(funx,(-self.slit.value/2),(self.slit.value/2))
		percent_l,percent_err_l = integrate.quad(funx,(-self.seeing.value/2),(self.seeing.value/2))
		percent = percent_u * percent_l # can use error if you add it later...
		extension = self.seeing.value * self.slit.value
		return percent,extension

	def recalculate_sky_flux(self,delta_lambda):
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
		if caller in etkeys.func_dichro:
			fblue_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y1, kind='cubic')
			fred_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y2, kind='cubic')
			self.stored_dichro_blue = fblue_dichro(self.wavelength_array)
			self.stored_dichro_red = fred_dichro(self.wavelength_array)
		return self.stored_dichro_blue,self.stored_dichro_red


	def recalculate_grating(self,caller):
		if caller in etkeys.func_grating:
			self.stored_grating_blue = spectres(self.wavelength_array,(dh.grating1[0]*10),dh.grating1[1])
			self.stored_grating_red = spectres(self.wavelength_array,(dh.grating2[0]*10),dh.grating2[1])
		return self.stored_grating_blue,self.stored_grating_red

	def recalculate_ccd(self,caller):
		if caller in etkeys.func_ccd:
			fblue_ccd = interpolate.interp1d((dh.ccd1[0]*10),dh.ccd1[1], kind='cubic')
			fred_ccd = interpolate.interp1d((dh.ccd2[0]*10),dh.ccd2[1], kind='cubic')
			self.stored_ccd_blue,self.stored_ccd_red = fblue_ccd(self.wavelength_array),fred_ccd(self.wavelength_array)
		return self.stored_ccd_blue,self.stored_ccd_red

	def recalculate_mirror_loss(self,caller):
		if caller in etkeys.func_mirror_loss:
			fmirror = interpolate.interp1d(dh.mirror_file_x,dh.mirror_file_y, kind='cubic')
			self.stored_mirror = fmirror(self.wavelength_array)
		return self.stored_mirror

	def recalculate_readnoise(self,caller): # probably not implemented in all necessary places yet...
		if caller in etkeys.func_readnoise:
			spectral_resolution = math.ceil((self.slit.value/(0.7/12))/2)*2 # px (ceil()/2)*2 to round up to next even integer
			spatial_resolution = math.ceil((self.seeing.value/(0.7/12))/2)*2
			extent = self.seeing.value * self.slit.value
			npix = extent/(0.7/12)**2
			self.stored_read_noise = math.ceil(dfs.rn_default * spectral_resolution * spatial_resolution / ((self.binning.active + 1)**2)) # +1 because zero indexing
			if not self.sss:
				print('{0} Pixel binning: ({1}x{1})'.format(dfs.string_prefix,(self.binning.active + 1)))
				print('{0} Extent: {1} arcsec^2\n{0} num pixels/resel: {2} px\n{0} spectral resolution: {3} px\n{0} spatial resolution: {4} px'.format(dfs.string_prefix,extent,int(math.ceil(npix)),spectral_resolution,spatial_resolution))
		return self.stored_read_noise

	def recalculate_atmospheric_extinction(self,caller):
		if caller in etkeys.func_atmo_ext:
			self.stored_extinction = spectres(self.wavelength_array,dh.atmo_ext_x,dh.atmo_ext_y)
		return self.stored_extinction