import os
import itertools
import numpy as np
import pandas as pd
import values as edl
import matplotlib.pyplot as plt
from spectres import spectres
from astropy import units as u
from astropy import constants as const
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.convolution import convolve
from scipy.integrate import quad

#obj,wavelength,filter_opt,magnitude,mag_sys_opt,grating_opt='LOW',redshift,exposure_time,seeing,slit_size,moon_day=0,plot_channel='BOTH',telescope_mode='FULL',binx=2,sss=False
class snr(object):
	def __init__(self,**kwargs): # wavelength,obj,filter_opt,slit_size,moon_day,mag_sys_opt,seeing,redshift
		try:
			self.__dict__ = dict(kwargs)
		except:
			print('NOPE.')

	def session(self):
		# required args
		self.seeing = self.seeing * u.arcsec
		self.slit_size = self.slit_size * u.arcsec
		self.wavelength = np.dot(self.wavelength,u.angstrom)
		
		# squelch
		try:
			sss = self.sss
		except:
			sss = False

		# object file
		try:
			# get custom pandas dataframe
			if isinstance(self.obj,pd.DataFrame):
				obj_x = self.obj[self.obj.columns[0]] # defaults to x = column 0
				obj_y = self.obj[self.obj.columns[1]] # defaults to y = column 1
			# or collect requested data file
			elif isinstance(self.obj,int):
				obj_panda = pd.read_csv(os.path.join(edl.galaxy_path,edl.galaxy_files[self.obj]),sep='\s+',skiprows=1)
				obj_x = obj_panda[obj_panda.columns[0]]
				obj_y = obj_panda[obj_panda.columns[1]]
			elif isinstance(self.obj,str):
				obj_panda = pd.read_csv(os.path.join(edl.galaxy_path,edl.galaxy_files[edl.galaxy_files.index(self.obj)]),sep='\s+',skiprows=1)
				obj_x = obj_panda[obj_panda.columns[0]]
				obj_y = obj_panda[obj_panda.columns[1]]
		except:
			print('[ error ] : You must supply an object (obj), either from our enumerated options, or a custom pandas dataframe.')

		# filter file
		try:
			# get custom pandas dataframe
			if isinstance(self.filter_opt,pd.DataFrame):
				sep_lambda = self.filter_opt[self.filter_opt.columns[0]]
				lambda_min = min(sep_lambda)
				lambda_max = max(sep_lambda)
			# or use one of our filter files
			elif isinstance(self.filter_opt,int):
				filter_panda = pd.read_csv(os.path.join(edl.filter_path,edl.filter_files[self.filter_opt]),sep=',',skiprows=1)
				filter_x = filter_panda[filter_panda.columns[0]]
				filter_y = filter_panda[filter_panda.columns[1]]
				lambda_min = min(filter_x)
				lambda_max = max(filter_x)
			elif isinstance(self.filter_opt,str):
				filter_panda = pd.read_csv(os.path.join(edl.filter_path,edl.filter_files[edl.filter_files.index(self.filter_opt)]),sep=',',skiprows=1)
				filter_x = filter_panda[filter_panda.columns[0]]
				filter_y = filter_panda[filter_panda.columns[1]]
				lambda_min = min(filter_x)
				lambda_max = max(filter_x)
		except:
			print('[ error ] : You must supply an filter (filter_opt), either from our enumerated options, or a custom pandas dataframe.')

		# skyfile
		try:
			if (int(self.moon_day) == 0):
				skyfile = pd.read_csv(os.path.join(edl.skyfiles_path,edl.skyfiles[0]),sep='\s+',skiprows=1)
			elif (int(self.moon_day) == 3):
				skyfile = pd.read_csv(os.path.join(edl.skyfiles_path,edl.skyfiles[1]),sep='\s+',skiprows=1)
			elif (int(self.moon_day) == 7):
				skyfile = pd.read_csv(os.path.join(edl.skyfiles_path,edl.skyfiles[2]),sep='\s+',skiprows=1)
			elif (int(self.moon_day) == 10):
				skyfile = pd.read_csv(os.path.join(edl.skyfiles_path,edl.skyfiles[3]),sep='\s+',skiprows=1)
			elif (int(self.moon_day) == 14):
				skyfile = pd.read_csv(os.path.join(edl.skyfiles_path,edl.skyfiles[4]),sep='\s+',skiprows=1)
		except:
			print('[ error ] : You must select an enumerated integer option \'number of days to new moon\' (moon_day)\n\tOptions: 0, 3, 7, 10, 14')
		# sky background
		sky_x = np.dot(skyfile[skyfile.columns[0]],(u.micron).to(u.angstrom))
		sky_y = np.dot(skyfile[skyfile.columns[1]],(1 / u.micron).to(1 / u.angstrom))
		old_res = sky_x[1] - sky_x[0]

		# area, permits string OR num mirrors
		try:
			if (self.telescope_mode.lower() == 'first') or (self.telescope_mode == 4): 
				area = edl.area[0] * (u.meter**2)
			elif (self.telescope_mode.lower() == 'full') or (self.telescope_mode == 7):
				area = edl.area[1] * (u.meter**2)
		except:
			area = edl.area[1] * (u.meter**2) # defaults to all theoretical mirrors
		
		# grating opt
		try:
			if (self.grating_opt.lower() == 'low'):
				delta_lambda_default = edl.dld[0] * (u.meter**2)
			elif (self.grating_opt.lower() == 'high'):
				delta_lambda_default = edl.dld[1] * (u.meter**2)
		except:
			delta_lambda_default = edl.dld[0] * (u.meter**2) # defaults to low resolution

		# atmospheric extinction
		atmo_ext = pd.read_csv(edl.atmo_ext_path,sep=',')
		self.end_extinction = (spectres(self.wavelength,np.asarray(atmo_ext[atmo_ext.columns[0]]),np.asarray(atmo_ext[atmo_ext.columns[1]])))
		
		# a couple other things
		lambda_a = np.dot(np.arange(lambda_min,lambda_max,0.1),u.angstrom)
		obj_x = obj_x * (1 + self.redshift) # apply redshift
		self.delta_lambda = delta_lambda_default * float(self.slit_size)/0.7
		self.plot_step = self.wavelength[-1]-self.wavelength[-2]
		
		if not sss:
			print('[ info ] : Delta lambda: {} Angstrom\nBinned pixel scale: {} Angstrom/pixel'.format(self.delta_lambda,self.plot_step))

		# calculate things
		self.flux_a = np.dot(spectres(lambda_a,obj_x,obj_y),u.jansky)
		self.trans = spectres(lambda_a,filter_x,filter_y)
		self.extinction = spectres(lambda_a,atmo_ext[atmo_ext.columns[0]],atmo_ext[atmo_ext.columns[1]])
		self.flux = [flux_a[i].cgs for i in range(flux_a)] # convert to cgs
		self._lambda = [lambda_a[i].to(u.meter) for i in range(lambda_a)] # convert to meter

		# validate flux
		blanks = itertools.count(0)
		try:
			for i in self.flux:
				if (self.flux[i] == 0):
					next(blanks)
				elif math.isnan(self.flux[i].value):
					self.flux[i] = 0
					next(blanks)
			if (blanks > (len(self.flux)/5)):
				zero_per = blanks / len(self.flux) * 100 # percent of values that are zero
				if not sss:
					print('[ error ] : Due to template limitations, {}% of spectrum in this bandpass is zero. Change redshift and/or bandpass.'.format(zero_per))
		except:
			if not sss:
				print('[ error ] : Due to template limitations, flux in this bandpass is zero. Change redshift and/or bandpass.')
			
		# magnitude system, vega or ab (default)
		try:
			if (self.mag_sys_opt.lower() == 'vega'):
				vega_panda = pd.read_csv(edl.vega_file,sep='\s+',skiprows=1) * u.jansky
				flux_vega = (spectres(lambda_a,vega_panda[vega_panda.columns[0]],vega_panda[vega_panda.columns[1]])).cgs
				self.mag_model = -2.5 * math.log(math.fsum((self.flux*self.trans*self.extinction*self._lambda)),10)/math.fsum(flux_vega*self.trans*self.extinction*self._lambda) + 0.026 # assuming that vega is 0.026 mag in all bands
			elif (self.mag_sys_opt.lower() == 'ab'):
				self.mag_model = -48.60 - 2.5 * math.log((math.fsum((self.flux*self.trans*self.extinction*self._lambda))/math.fsum((self.trans*self._lambda*self.extinction*(const.c/(self._lambda**2))))),10) # zeropoint is 48.60
		except:
			print('[ error ] : Attribute \'mag_sys_opt\' must be set! Include \'VEGA\' or \'AB\'')

		self.del_mag = self.magnitude - self.mag_model # desired vs template
		self.ss_lambda = obj_x # for clarity
		self.ss_flux = obj_y * 10 ** (-self.del_mag/2.5) # rescale

		old_res = self.ss_lambda[2] - self.ss_lambda[1]
		if (old_res < self.plot_step):
			sigma = self.delta_lambda / gaussian_sigma_to_fwhm # delta_lambda is a scalar too
			_x = np.arange(start=(-5*sigma),stop=(5*sigma),step=old_res) # therefore, this
			funx = (1/(sigma*(math.sqrt(2*math.pi))))*math.exp(-_x**2/(2*sigma**2))
			kernel = funx / np.trapz(funx)
			self.ss_flux = convolve(self.ss_flux,kernel,'fill')
		else:
			print('[ error ] : The object template resolution is lower than the instrument\'s resolution.')

		ts_flux = (spectres(self.wavelength,self.ss_lambda,self.ss_flux)) * u.jansky
		power = ((ts_flux * 1e-3) * area * self.exposure_time * self.plot_step) # should come out as joules
		counts_counts = (power / ((const.h * const.c) / self.wavelength)) / 1e10
		sigma = self.seeing / gaussian_sigma_to_fwhm
		funx = lambda _x : (1/(sigma*(math.sqrt(2 * math.pi))) * math.exp(-_x**2 / (2 * sigma**2)))
		self.percent = quad(funx,(-self.slit_size/2),(self.slit_size/2)) * quad(funx,(-self.seeing/2),(self.seeing/2)) # for signal calculation later
		self.extension = (self.seeing * self.slit_size)

		# calculate noise
		sigma = self.delta_lambda / gaussian_sigma_to_fwhm
		_x = np.arange(start=(-5*sigma),stop=(5*sigma),step=old_res)
		funx = (1/(sigma*(math.sqrt(2*math.pi))))*math.exp(-_x**2/(2*sigma**2))
		kernel = funx / np.trapz(funx)
		self.sky_flux = convolve(sky_y,kernel,'fill') * u.jansky
		self.counts_noise = (self.sky_flux * extension * area) # number of photons. should be unitless
		
		# resample dichroic transmission
		dichroic_panda = pd.read_csv(edl.dichroic_path,sep='\s+',skiprows=1)
		dichro_x = (dichroic_panda[dichroic_panda.columns[0]])
		dichro_x = [(dichro_x[i] * u.nm).to(u.angstrom) for i in range (dichro_x)]
		dichro_y1 = (dichroic_panda[dichroic_panda.columns[1]]) # reflectivity ... units? joule per square meter or something?
		dichro_y2 = (dichroic_panda[dichroic_panda.columns[2]]) # transmission
		self.red_dichro = (spectres(self.wavelength,dichro_x,dichro_y1))
		self.blue_dichro = (spectres(self.wavelength,dichro_x,dichro_y2))

		# resample grating efficiency
		grating_panda_red = pd.read_csv(edl.grating_path[0],sep='\s+',skiprows=1)
		grating_panda_blue = pd.read_csv(edl.grating_path[1],sep='\s+',skiprows=1)
		self.red_grating = (spectres(self.wavelength,grating_panda_red[grating_panda_red.columns[0]],grating_panda_red[grating_panda_red.columns[1]]))
		self.blue_grating = (spectres(self.wavelength,grating_panda_blue[grating_panda_blue.columns[0]],grating_panda_blue[grating_panda_blue.columns[1]]))

		# resample ccd efficiency
		ccd_panda_red = pd.read_csv(edl.ccd_path[0],sep='\s+',skiprows=1)
		ccd_panda_blue = pd.read_csv(edl.ccd_path[1],sep='\s+',skiprows=1)
		self.red_ccd = (spectres(self.wavelength,ccd_panda_red[ccd_panda_red.columns[0]],ccd_panda_red[ccd_panda_red.columns[1]]))
		self.blue_ccd = (spectres(self.wavelength,ccd_panda_blue[ccd_panda_blue.columns[0]],ccd_panda_blue[ccd_panda_blue.columns[1]]))

		# read noise
		self.spectral_resolution = math.ceil((self.slit_size/(0.7/12))/2)*2 #px (ceil()/2)*2 to round up to next even integer
		self.spatial_resolution = math.ceil((self.seeing/(0.7/12))/2)*2 #px (ceil()/2)*2 to round up to next even integer 
		self.extent = self.seeing * self.slit_size
		self.npix = self.extent/(0.7/12)**2
		try:
			self.rn = float(self.rn)
		except:
			self.rn = edl.rn_default
		try:
			if isinstance(self.bin_size,int):
				if (self.bin_size > 0) and (self.bin_size < 5):
					if not sss:
						pass
			elif (self.bin_size in edl.bin_options_str):
				self.bin_size = int(self.bin_size[0])
			elif (self.bin_size == 'default'):
				self.bin_size = edl.bin_options_int[edl.bin_options_default_index]
		except:
			self.bin_size = edl.bin_options_int[edl.bin_options_default_index]
		
		print('[ etc ] : Binning ({}x{})'.format(self.bin_size,self.bin_size))
		self.readnoise = int(self.rn * self.spectral_resolution * self.spatial_resolution / (self.bin_size**2))
		if not sss:
			print('extent: {} arcsec^2\nnum pixels: {} px\nspectral resolution: {} px\nspatial resolution: {} px'.format(self.extent,self.npix,self.spectral_resolution,self.spatial_resolution))

		''' calculations '''
		# signal
		self.red_total_efficiency = (self.red_dichro * self.red_grating * self.red_ccd * edl.coating_efficiency * self.end_extinction)
		self.blue_total_efficiency = (self.blue_dichro * self.blue_grating * self.blue_ccd * edl.coating_efficiency * self.end_extinction)
		self.red_signal = (self.counts_counts * self.percent * self.red_total_efficiency)
		self.blue_signal = (self.counts_counts * self.percent * self.blue_total_efficiency)

		# noise
		self.red_total_efficiency_noise = (self.red_dichro * self.red_grating * self.red_ccd * edl.coating_efficiency)
		self.blue_total_efficiency_noise = (self.blue_dichro * self.blue_grating * self.blue_ccd * edl.coating_efficiency)
		self.red_noise = (self.counts_noise * self.red_total_efficiency_noise)
		self.blue_noise = (self.counts_noise * self.blue_total_efficiency_noise)

		# signal/noise ratio
		self.red_snr = (self.red_signal / np.sqrt(red_signal + red_noise + (self.readnoise**2))) # add thing here for ummm binning
		self.blue_snr = (self.blue_signal / np.sqrt(blue_signal + blue_noise + (self.readnoise**2)))

		# generate error based on signal/noise ratio
		self.e_norm_r = np.random.normal(loc=0, scale=(np.sqrt(self.red_signal + self.red_noise)),size=len(self.red_snr))
		self.e_norm_b = np.random.normal(loc=0, scale=(np.sqrt(self.blue_signal + self.blue_noise)),size=len(self.blue_snr))

		# plot-ready array objects
		self.snr = {self.red_snr, self.blue_snr}
		self.obs_spec = {self.red_signal, self.blue_signal}
		self.obs_specnoise = {np.add(self.red_signal,self.e_norm_r),np.add(self.blue_signal,self.e_norm_b)}
		self.osb_sky_background = {self.red_noise,self.blue_noise}

	def try_plotting(self):
		plt.plot(self.wavelength,self.snr)
		plt.show()