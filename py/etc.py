import os
import math
import itertools
import numpy as np
import pandas as pd
import values as edl
from spectres import spectres
from astropy import units as u
from astropy import constants as const
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.convolution import convolve
from scipy.integrate import quad

#obj,wavelength,filter_opt,magnitude,mag_sys_opt,grating_opt='LOW',redshift,exposure_time,seeing,slit_size,moon_days=0,plot_channel='BOTH',telescope_mode='FULL',binx=2,sss=False
class simulate:
	def __init__(self,obj,filter_opt,wavelength,seeing,slit_size,moon_days,redshift,**kwargs): # wavelength,obj,filter_opt,slit_size,moon_days,mag_sys_opt,seeing,redshift
		self.__dict__ = dict(kwargs) # pick up any optionals
		# required args
		self.redshift = redshift
		self.obj = obj
		self.filter_opt = filter_opt
		#self.wavelength = np.dot(wavelength,u.angstrom)
		self.wavelength = wavelength
		self.seeing = seeing # * u.arcsec
		self.slit_size = slit_size # * u.arcsec
		self.moon_days = moon_days
		self.delta_lambda_default = edl.dld[0] # * (u.meter**2) # default low res
		# subsiste sermonem statim
		try:
			sss = self.sss
		except:
			self.sss = False



	# grating opt
	def grating_opt(self,_grating_opt=None):
		try:
			if (_grating_opt == None):
				_grating_opt = self.grating_opt # public to private
			else:
				pass # use public as private
		except:
			if (_grating_opt != None):
				if (_grating_opt.lower() == 'low'):
					self.grating_opt = _grating_opt # update public
					self.delta_lambda_default = edl.dld[0] # * (u.meter**2)
				elif (_grating_opt.lower() == 'high'):
					self.grating_opt = _grating_opt # update public
					self.delta_lambda_default = edl.dld[1] # * (u.meter**2)
				else:
					print('[ error ] : You must specify a valid grating mode!\nOptions are: \'high\' and \'low\' (default)')
			else:
				pass			


	# get object file
	def objectfile(self):
		try:
			# collect requested data file
			if isinstance(self.obj,int):
				obj_panda = pd.read_csv(os.path.join(edl.galaxy_path,edl.galaxy_files[self.obj]),sep='\s+',skiprows=1)
				obj_x = obj_panda[obj_panda.columns[0]]
				obj_y = obj_panda[obj_panda.columns[1]]
			elif isinstance(self.obj,str):
				obj_panda = pd.read_csv(os.path.join(edl.galaxy_path,edl.galaxy_files[edl.galaxy_files.index(self.obj)]),sep='\s+',skiprows=1)
				obj_x = obj_panda[obj_panda.columns[0]]
				obj_y = obj_panda[obj_panda.columns[1]]
			# get custom pandas dataframe
			elif isinstance(self.obj,pandas.core.frame.DataFrame):
				obj_x = self.obj[self.obj.columns[0]] # defaults to x = column 0
				obj_y = self.obj[self.obj.columns[1]] # defaults to y = column 1

			self.obj_x = obj_x
			self.obj_y = obj_y
		except:
			print('[ error ] : You must supply an object (obj), either from our enumerated options, or a custom pandas dataframe.')


	# get filter file
	def filterfile(self):
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

		#self.lambda_a = np.dot(np.arange(lambda_min,lambda_max,0.1),u.angstrom)
		self.lambda_a = np.arange(lambda_min,lambda_max,0.1)
		obj_x = self.obj_x * (1 + self.redshift) # apply redshift
		self.delta_lambda = self.delta_lambda_default * self.slit_size/0.7
		self.plot_step = self.wavelength[-1]-self.wavelength[-2]
		if not self.sss:
			print('[ info ] : Delta lambda: {} Angstrom\n[ info ] : Binned pixel scale: {} Angstrom/pixel'.format(self.delta_lambda,self.plot_step))


	# atmospheric extinction
	def atmoext(self):
		try:
			self.end_extinction = spectres(self.wavelength,np.asarray(atmo_ext[atmo_ext.columns[0]]),np.asarray(atmo_ext[atmo_ext.columns[1]]))
		except:
			atmo_ext = pd.read_csv(edl.atmo_ext_path,sep=',')
			self.end_extinction = spectres(self.wavelength,np.asarray(atmo_ext[atmo_ext.columns[0]]),np.asarray(atmo_ext[atmo_ext.columns[1]]))


	# skyfile depends on number of days since/to a new moon
	def skyfile(self,_moon_day=None):
		try:
			if (_moon_day == None):
				_moon_day = self.moon_days # public to private
			else:
				pass # use public as private
		except:
			if ( _moon_day != None):
				self.moon_days = _moon_day # private to public
			else:
				pass # use private as private

		if (int(_moon_day) == 0):
			skyfile = pd.read_csv(os.path.join(edl.skyfiles_path,edl.skyfiles[0]),sep='\s+',skiprows=1)
		elif (int(_moon_day) == 3):
			skyfile = pd.read_csv(os.path.join(edl.skyfiles_path,edl.skyfiles[1]),sep='\s+',skiprows=1)
		elif (int(_moon_day) == 7):
			skyfile = pd.read_csv(os.path.join(edl.skyfiles_path,edl.skyfiles[2]),sep='\s+',skiprows=1)
		elif (int(_moon_day) == 10):
			skyfile = pd.read_csv(os.path.join(edl.skyfiles_path,edl.skyfiles[3]),sep='\s+',skiprows=1)
		elif (int(_moon_day) == 14):
			skyfile = pd.read_csv(os.path.join(edl.skyfiles_path,edl.skyfiles[4]),sep='\s+',skiprows=1)
		# sky background unit conversions
		sky_x = (np.dot(skyfile[skyfile.columns[0]],(u.micron).to(u.angstrom)))
		self.sky_y = (np.dot(skyfile[skyfile.columns[1]],(1 / u.micron).to(1 / u.angstrom)))
		self.old_res = sky_x[1] - sky_x[0]


	# area, permits string OR num mirrors
	def mirrors(self,_telescope_mode=None):
		try:
			if (_telescope_mode == None):
				_telescope_mode = self.telescope_mode # public to private
			else:
				pass # use public as private
		except:
			if (_telescope_mode != None):
				self.telescope_mode = _telescope_mode # private to public
			else:
				pass # use private as private

		try:
			#if (_telescope_mode.lower() == 'first') or (_telescope_mode == 4): gotta check type first
			#elif (_telescope_mode.lower() == 'full') or (_telescope_mode == 7):
			if (_telescope_mode == 4): 
				self.area = edl.area[0] # * (u.meter**2)
			elif (_telescope_mode == 7):
				self.area = edl.area[1] # * (u.meter**2)
		except:
			print('[ error ] : You must specify a valid number of mirrors (telescope_mode)\nOptions are: 4 (first light) and 7 (full size)')


	# resample dichroic transmission
	def resample_dichroic(self,new_wavelength=None):
		try:
			if (new_wavelength == None):
				new_wavelength = self.wavelength # use public
			else:
				pass
		except:
			if (new_wavelength != None): # use private
				if (isinstance(new_wavelength,numpy.ndarray)):
					if (len(new_wavelength) <= 20): # probably a mistake
						print('[ error ] : Resampling dichroic transmission, new wavelength appears to be invalid. Perhaps you meant to specify start/end/step in a tuple?')
					else:
						self.wavelength = new_wavelength # update public
				elif (isinstance(new_wavelength,tuple)):
					if (len(new_wavelength) == 3):
						new_wavelength = np.arange(new_wavelength[0],new_wavelength[1],new_wavelength[2]) # build new wavelength one from tuple
						self.wavelength = new_wavelength # update public
					else: # uh oh
						print('[ error ] : Resampling dichroic transmission, tuple wavelength must be within valid ranges, and in format (start,stop,step).')

		try: # problem here :(
			#self.wavelength = new_wavelength
			self.red_dichro = spectres(new_wavelength,dichro_x,dichro_y1)
			self.blue_dichro = spectres(new_wavelength,dichro_x,dichro_y2)
		except:
			dichroic_panda = pd.read_csv(edl.dichroic_path,sep='\s+',skiprows=1)
			dichro_x = np.asarray(dichroic_panda[dichroic_panda.columns[0]])
			#dichro_x = [(dichro_x[i] * u.nm).to(u.angstrom) for i in range(len(dichro_x))]
			dichro_x = np.dot(dichro_x,10) # until I learn astropy a bit better...
			dichro_y1 = np.asarray(dichroic_panda[dichroic_panda.columns[1]]) # reflectivity ... units? joule per square meter or something?
			dichro_y2 = np.asarray(dichroic_panda[dichroic_panda.columns[2]]) # transmission
			#self.red_dichro = spectres(new_wavelength,dichro_x,dichro_y1)
			#self.blue_dichro = spectres(new_wavelength,dichro_x,dichro_y2)
			#self.red_dichro = convolve(,new_wavelength)


	# resample grating efficiency
	def resample_grating(self,new_wavelength=None):
		new_wavelength = self.wavelength # debugging...
		grating_panda_red = pd.read_csv(edl.grating_path[0],sep=',',skiprows=1)
		grating_panda_blue = pd.read_csv(edl.grating_path[1],sep=',',skiprows=1)
		grating_red_x = np.dot(np.asarray(grating_panda_red[grating_panda_red.columns[0]]),10)
		grating_red_y = np.dot(np.asarray(grating_panda_red[grating_panda_red.columns[1]]),10)
		grating_blue_x = np.dot(np.asarray(grating_panda_blue[grating_panda_blue.columns[0]]),10)
		grating_blue_y = np.dot(np.asarray(grating_panda_blue[grating_panda_blue.columns[1]]),10)
		self.red_grating = spectres(new_wavelength,grating_red_x,grating_red_y)
		self.blue_grating = spectres(new_wavelength,grating_blue_x,grating_blue_y)


	# resample ccd efficiency
	def resample_ccd(self):
		new_wavelength = self.wavelength # debugging...
		try:
			self.red_ccd = spectres(new_wavelength,ccd_panda_red[ccd_panda_red.columns[0]],ccd_panda_red[ccd_panda_red.columns[1]])
			self.blue_ccd = spectres(new_wavelength,ccd_panda_blue[ccd_panda_blue.columns[0]],ccd_panda_blue[ccd_panda_blue.columns[3]])
		except:
			ccd_panda_red = pd.read_csv(edl.ccd_path[0],sep='\s+',skiprows=1)
			ccd_panda_blue = pd.read_csv(edl.ccd_path[1],sep='\s+',skiprows=1)
			ccd_red_x = np.dot(np.asarray(ccd_panda_red[ccd_panda_red.columns[0]]),10)
			ccd_red_y = np.dot(np.asarray(ccd_panda_red[ccd_panda_red.columns[1]]),10)
			ccd_blue_x = np.dot(np.asarray(ccd_panda_blue[ccd_panda_blue.columns[0]]),10)
			ccd_blue_y = np.dot(np.asarray(ccd_panda_blue[ccd_panda_blue.columns[3]]),10)
			self.red_ccd = spectres(new_wavelength,ccd_red_x,ccd_red_y)
			self.blue_ccd = spectres(new_wavelength,ccd_blue_x,ccd_blue_y)


	def readnoise(self):
		# read noise
		self.spectral_resolution = math.ceil((self.slit_size/(0.7/12))/2)*2 #px (ceil()/2)*2 to round up to next even integer
		self.spatial_resolution = math.ceil((self.seeing/(0.7/12))/2)*2 #px (ceil()/2)*2 to round up to next even integer 
		self.extent = self.seeing * self.slit_size
		self.npix = self.extent/(0.7/12)**2
		# default unless specified
		try:
			self.rn = float(self.rn)
		except:
			self.rn = edl.rn_default
		try:
			if isinstance(self.bin_size,int):
				if (self.bin_size > 0) and (self.bin_size < 5):
					if not self.sss:
						pass
			elif (self.bin_size in edl.bin_options_str):
				self.bin_size = int(self.bin_size[0])
			elif (self.bin_size == 'default'):
				self.bin_size = edl.bin_options_int[edl.bin_options_default_index]
		except:
			self.bin_size = edl.bin_options_int[edl.bin_options_default_index]
		
		self.readnoise = int(self.rn * self.spectral_resolution * self.spatial_resolution / (self.bin_size**2))
		if not self.sss:
			print('[ info ] : Binning ({}x{})'.format(self.bin_size,self.bin_size))
			print('[ info ] : Extent: {} arcsec^2\n[ info ] : num pixels: {} px\n[ info ] : spectral resolution: {} px\n[ info ] : spatial resolution: {} px'.format(self.extent,self.npix,self.spectral_resolution,self.spatial_resolution))


	def resample_flux(self):
		# resample flux
		#self.flux_a = np.dot(spectres(self.lambda_a,obj_x,obj_y),u.jansky)
		self.flux_a = spectres(self.lambda_a,obj_x,obj_y)
		self.trans = spectres(self.lambda_a,filter_x,filter_y)
		self.extinction = spectres(self.lambda_a,atmo_ext[atmo_ext.columns[0]],atmo_ext[atmo_ext.columns[1]])
		#self.flux = [flux_a[i].cgs for i in range(flux_a)] # convert to cgs
		self.flux = np.dot(self.flux_a,1e10)
		#self._lambda = [self.lambda_a[i].to(u.meter) for i in range(self.lambda_a)] # convert to meter
		self._lambda = np.divide(self.lambda_a,1e10)
		# validate resampled flux
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
				print('[ error ] : Due to template limitations, {}% of spectrum in this bandpass is zero. Change redshift and/or bandpass.'.format(zero_per))
		except:
			print('[ error ] : Due to template limitations, flux in this bandpass is zero. Change redshift and/or bandpass.')
		# compare resampled flux to object's template
		self.del_mag = self.magnitude - self.mag_model # desired vs template
		self.ss_lambda = obj_x
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


	# magnitude system, vega or ab (default)
	def magsys(self):
		#try:
		if (self.mag_sys_opt.lower() == 'vega'):
			vega_panda = pd.read_csv(edl.vega_file,sep='\s+',skiprows=1) # * u.jansky
			flux_vega = np.dot(spectres(self.lambda_a,np.asarray(vega_panda[vega_panda.columns[0]]),np.asarray(vega_panda[vega_panda.columns[1]])),1e10) #.cgs
			self.mag_model = -2.5 * math.log(math.fsum((self.flux*self.trans*self.extinction*self._lambda)),10)/math.fsum(flux_vega*self.trans*self.extinction*self._lambda) + 0.026 # assuming that vega is 0.026 mag in all bands
		elif (self.mag_sys_opt.lower() == 'ab'):
			self.mag_model = -48.60 - 2.5 * math.log((math.fsum((self.flux*self.trans*self.extinction*self._lambda))/math.fsum((self.trans*self._lambda*self.extinction*(const.c/(self._lambda**2))))),10) # zeropoint is 48.60
		#except:
		#	print('[ error ] : Attribute \'mag_sys_opt\' must be set! Include \'VEGA\' or \'AB\'')


	def signal(self):
		# counts counts
		ts_flux = (spectres(self.wavelength,self.ss_lambda,self.ss_flux)) # * u.jansky
		power = ((ts_flux * 1e-3) * self.area * self.exposure_time * self.plot_step) # should come out as joules
		counts_counts = (power / ((const.h * const.c) / self.wavelength)) / 1e10
		sigma = self.seeing / gaussian_sigma_to_fwhm
		funx = lambda _x : (1/(sigma*(math.sqrt(2 * math.pi))) * math.exp(-_x**2 / (2 * sigma**2)))
		self.percent = quad(funx,(-self.slit_size/2),(self.slit_size/2)) * quad(funx,(-self.seeing/2),(self.seeing/2)) # for signal calculation later
		self.extension = (self.seeing * self.slit_size)
		# final step
		self.red_total_efficiency = (self.red_dichro * self.red_grating * self.red_ccd * edl.coating_efficiency * self.end_extinction)
		self.blue_total_efficiency = (self.blue_dichro * self.blue_grating * self.blue_ccd * edl.coating_efficiency * self.end_extinction)
		self.red_signal = (self.counts_counts * self.percent * self.red_total_efficiency)
		self.blue_signal = (self.counts_counts * self.percent * self.blue_total_efficiency)


	def noise(self):
		# calculate noise
		sigma = self.delta_lambda / gaussian_sigma_to_fwhm
		_x = np.arange(start=(-5*sigma),stop=(5*sigma),step=old_res)
		funx = (1/(sigma*(math.sqrt(2*math.pi))))*math.exp(-_x**2/(2*sigma**2))
		kernel = funx / np.trapz(funx)
		self.sky_flux = convolve(sky_y,kernel,'fill') # * u.jansky
		self.counts_noise = (self.sky_flux * extension * self.area) # number of photons
		# final step
		self.red_total_efficiency_noise = (self.red_dichro * self.red_grating * self.red_ccd * edl.coating_efficiency)
		self.blue_total_efficiency_noise = (self.blue_dichro * self.blue_grating * self.blue_ccd * edl.coating_efficiency)
		self.red_noise = (self.counts_noise * self.red_total_efficiency_noise)
		self.blue_noise = (self.counts_noise * self.blue_total_efficiency_noise)


	# signal/noise ratio
	def calc_snr(self):
		self.red_snr = (self.red_signal / np.sqrt(red_signal + red_noise + (self.readnoise**2)))
		self.blue_snr = (self.blue_signal / np.sqrt(blue_signal + blue_noise + (self.readnoise**2)))


	# generate error based on signal/noise ratio
	def gen_err(self):
		self.e_norm_r = np.random.normal(loc=0, scale=(np.sqrt(self.red_signal + self.red_noise)),size=len(self.red_snr))
		self.e_norm_b = np.random.normal(loc=0, scale=(np.sqrt(self.blue_signal + self.blue_noise)),size=len(self.blue_snr))

#	def state(self,caller,command_index):
#		if (caller == 0): # snr
	def snr(self):
		self.grating_opt()
		self.objectfile()
		self.filterfile()
		self.skyfile()
		self.mirrors()
		self.resample_dichroic()
		self.resample_grating()
		self.resample_ccd()
		self.readnoise()
		self.resample_flux()
		self.magsys()
		self.signal()
		self.noise()
		self.calc_snr()
		self.gen_err()
		return self.red_snr, self.blue_snr


	def obs_spec(self,noise=True):
		if noise:
			return np.add(self.red_signal,self.e_norm_r),np.add(self.blue_signal,self.e_norm_b)
		else:
			return self.red_signal, self.blue_signal


	def sky_background(self):
		return self.red_noise,self.blue_noise


'''
class manifest:
	# to implement
'''

'''
class specs:
	def dichroic_throughput():
		pass

	
	def grating_throughput():
		pass


	def ccd_qe():
		pass


	def atmospheric_extinction():
		try:
			atmo_ext = spectres(self.wavelength,np.asarray(atmo_ext[atmo_ext.columns[0]]),np.asarray(atmo_ext[atmo_ext.columns[1]]))
		except:
			atmo_ext_panda = pd.read_csv(edl.atmo_ext_path,sep=',')
			atmo_ext = spectres(self.wavelength,np.asarray(atmo_ext[atmo_ext.columns[0]]),np.asarray(atmo_ext[atmo_ext.columns[1]]))

'''