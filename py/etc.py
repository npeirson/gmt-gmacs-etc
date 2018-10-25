import os
import math
import numpy as np
import pandas as pd
import values as edl
from spectres import spectres
from astropy import units as u
from astropy import constants as const
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel
from scipy.integrate import quad, quadrature


#obj,wavelength,filter_opt,magnitude,mag_sys_opt,grating_opt='LOW',redshift,exposure_time,seeing,slit_size,moon_days=0,plot_channel='BOTH',telescope_mode='FULL',binx=2,sss=False
class simulate:
	def __init__(self,obj,filter_opt,magnitude,wavelength,seeing,slit_size,moon_days,redshift,sss=False,telescope_mode=4,**kwargs): # wavelength,obj,filter_opt,slit_size,moon_days,mag_sys_opt,seeing,redshift
		self.__dict__ = dict(kwargs) # pick up any optionals
		# required args
		self.magnitude = magnitude
		self.redshift = redshift
		self.obj = obj
		self.filter_opt = filter_opt
		self.wavelength = np.arange(8000,10000,10)
		self.telescope_mode = telescope_mode
		self.seeing = seeing # * u.arcsec
		self.slit_size = slit_size # * u.arcsec
		self.moon_days = moon_days
		self.delta_lambda_default = edl.dld[0] # * (u.meter**2) # default low res
		self.sss = sss	


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
			isinstance(self.object_type)
		except:
			self.object_type = 'stellar'

		if (self.object_type.lower() == 'galaxy') or (self.object_type.lower() == 'galactic'):
			folder_path = edl.galaxy_path
			file_path = edl.galaxy_files
			print('[ info ] : Object type: galactic')
		else:
			folder_path = edl.stellar_path
			file_path = edl.stellar_files
			print('[ info ] : Object type: stellar')

		try:
			# collect requested data file
			if isinstance(self.obj,int):
				obj_panda = pd.read_csv(os.path.join(folder_path,file_path[self.obj]),sep='\s+',skiprows=1)
				obj_x = obj_panda[obj_panda.columns[0]]
				obj_y = obj_panda[obj_panda.columns[1]]
				print('[ info ] : Model object: {}'.format(file_path[self.obj][:-4]))
			elif isinstance(self.obj,str):
				obj_panda = pd.read_csv(os.path.join(folder_path,file_path[file_path.index(self.obj)]),sep='\s+',skiprows=1)
				obj_x = obj_panda[obj_panda.columns[0]]
				obj_y = obj_panda[obj_panda.columns[1]]
			# get custom pandas dataframe
			elif isinstance(self.obj,pandas.core.frame.DataFrame):
				obj_x = self.obj[self.obj.columns[0]] # defaults to x = column 0
				obj_y = self.obj[self.obj.columns[1]] # defaults to y = column 1

			self.obj_x = np.dot(np.asarray(obj_x),(1 + self.redshift)) # apply redshift
			self.obj_y = np.asarray(obj_y)
		except:
			print('[ error ] : You must supply an object (obj), either from our enumerated options, or a custom pandas dataframe.')


	# get filter file
	def filterfile(self):
		try:
			if isinstance(self.filter_opt,pd.DataFrame): # get custom pandas dataframe
				sep_lambda = self.filter_opt[self.filter_opt.columns[0]]
				lambda_min = min(sep_lambda)
				lambda_max = max(sep_lambda)
			# or use one of our filter files
			elif isinstance(self.filter_opt,int):
				filter_panda = pd.read_csv(os.path.join(edl.filter_path,edl.filter_files[self.filter_opt]),sep=',')
				filter_x = filter_panda[filter_panda.columns[0]]
				filter_y = filter_panda[filter_panda.columns[1]]
				lambda_min = min(filter_x)
				lambda_max = max(filter_x)
			elif isinstance(self.filter_opt,str):
				filter_panda = pd.read_csv(os.path.join(edl.filter_path,edl.filter_files[edl.filter_files.index(self.filter_opt)]),sep=',')
				filter_x = filter_panda[filter_panda.columns[0]]
				filter_y = filter_panda[filter_panda.columns[1]]
				lambda_min = min(filter_x)
				lambda_max = max(filter_x)

			self.filter_x = filter_x
			self.filter_y = filter_y
			self.lambda_a = np.arange(lambda_min,lambda_max,0.1)
			self.delta_lambda = self.delta_lambda_default * self.slit_size/0.7
			self.plot_step = self.wavelength[-1] - self.wavelength[-2]
			if self.sss:
				print('[ info ] : Delta lambda: {} Angstrom\n[ info ] : Binned pixel scale: {} Angstrom/pixel'.format(self.delta_lambda,self.plot_step))
		
		except:
			print('[ error ] : You must supply an filter (filter_opt), either from our enumerated options, or a custom pandas dataframe.')


	# atmospheric extinction
	def atmoext(self):
		try:
			self.end_extinction = spectres(self.wavelength,np.asarray(self.atmo_ext[self.atmo_ext.columns[0]]),np.asarray(self.atmo_ext[self.atmo_ext.columns[1]]))
		except:
			self.atmo_ext = pd.read_csv(edl.atmo_ext_path,sep=',')
			self.end_extinction = spectres(self.wavelength,np.asarray(self.atmo_ext[self.atmo_ext.columns[0]]),np.asarray(self.atmo_ext[self.atmo_ext.columns[1]]))


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
		self.sky_x = np.dot(skyfile[skyfile.columns[0]],10) # micron to Angstrom
		self.sky_y = np.dot(skyfile[skyfile.columns[1]],10) # micron to Angstrom
		self.old_res_sky = self.sky_x[1] - self.sky_x[0]


	# area, permits string OR num mirrors
	def mirrors(self):
		_telescope_mode = self.telescope_mode # debugging
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
			if (_telescope_mode == 4): 
				self.area = edl.area[0] # * (u.meter**2)
			elif (_telescope_mode == 7):
				self.area = edl.area[1] # * (u.meter**2)
		except:
			print('[ error ] : You must specify a valid number of mirrors (telescope_mode)\nOptions are: 4 (first light) and 7 (full size)')


	# resample dichroic transmission
	def resample_dichroic(self):
		dichroic_panda = pd.read_csv(edl.dichroic_path,sep='\s+',skiprows=1)
		dichro_x = np.asarray(dichroic_panda[dichroic_panda.columns[0]])
		dichro_x = np.dot(dichro_x,10) # nm to Angstrom
		dichro_y1 = np.asarray(dichroic_panda[dichroic_panda.columns[1]]) # reflectivity ... units? joule per square meter or something?
		dichro_y2 = np.asarray(dichroic_panda[dichroic_panda.columns[2]]) # transmission
		self.red_dichro = spectres(self.wavelength,dichro_x,dichro_y1)
		self.blue_dichro = spectres(self.wavelength,dichro_x,dichro_y2)
			

	# resample grating efficiency
	def resample_grating(self):
		grating_panda_red = pd.read_csv(edl.grating_path[0],sep='\s+')
		grating_panda_blue = pd.read_csv(edl.grating_path[1],sep='\s+')
		grating_red_x = np.dot(np.asarray(grating_panda_red[grating_panda_red.columns[0]]),10)
		grating_red_y = np.dot(np.asarray(grating_panda_red[grating_panda_red.columns[1]]),10)
		grating_blue_x = np.dot(np.asarray(grating_panda_blue[grating_panda_blue.columns[0]]),10)
		grating_blue_y = np.dot(np.asarray(grating_panda_blue[grating_panda_blue.columns[1]]),10)
		if (self.plot_channel == 'red'):
			self.red_grating = spectres(self.wavelength,grating_red_x,grating_red_y)
			self.blue_grating = spectres(self.wavelength,grating_blue_x,grating_blue_y)


	# resample ccd efficiency
	def resample_ccd(self):
		try:
			self.red_ccd = spectres(self.wavelength,ccd_panda_red[ccd_panda_red.columns[0]],ccd_panda_red[ccd_panda_red.columns[1]])
			self.blue_ccd = spectres(self.wavelength,ccd_panda_blue[ccd_panda_blue.columns[0]],ccd_panda_blue[ccd_panda_blue.columns[3]])
		except:
			ccd_panda_red = pd.read_csv(edl.ccd_path[0],sep='\s+',skiprows=1)
			ccd_panda_blue = pd.read_csv(edl.ccd_path[1],sep='\s+',skiprows=1)
			ccd_red_x = np.dot(np.asarray(ccd_panda_red[ccd_panda_red.columns[0]]),10)
			ccd_red_y = np.dot(np.asarray(ccd_panda_red[ccd_panda_red.columns[1]]),10)
			ccd_blue_x = np.dot(np.asarray(ccd_panda_blue[ccd_panda_blue.columns[0]]),10)
			ccd_blue_y = np.dot(np.asarray(ccd_panda_blue[ccd_panda_blue.columns[3]]),10)
			self.red_ccd = spectres(self.wavelength,ccd_red_x,ccd_red_y)
			self.blue_ccd = spectres(self.wavelength,ccd_blue_x,ccd_blue_y)


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
						print('[ info ] : Pixel binning: ({}x{})'.format(self.bin_size,self.bin_size))
					else:
						pass
			elif (self.bin_size in edl.bin_options_str):
				self.bin_size = int(self.bin_size[0])
			elif (self.bin_size == 'default'):
				self.bin_size = edl.bin_options_int[edl.bin_options_default_index]
		except:
			self.bin_size = edl.bin_options_int[edl.bin_options_default_index]
		
		self.readnoise = math.ceil(self.rn * self.spectral_resolution * self.spatial_resolution / (self.bin_size**2))
		if not self.sss:
			print('[ info ] : Binning ({}x{})'.format(self.bin_size,self.bin_size))
			print('[ info ] : Extent: {} arcsec^2\n[ info ] : num pixels: {} px\n[ info ] : spectral resolution: {} px\n[ info ] : spatial resolution: {} px'.format(self.extent,self.npix,self.spectral_resolution,self.spatial_resolution))


	def resample_flux(self):
		# resample flux
		self.flux_a = spectres(self.wavelength,self.obj_x,self.obj_y)
		self.trans = np.trapz(np.asarray(self.filter_y),np.asarray(self.filter_x))
		self.extinction = spectres(self.wavelength,np.asarray(self.atmo_ext[self.atmo_ext.columns[0]]),np.asarray(self.atmo_ext[self.atmo_ext.columns[1]]))
		self.flux = np.dot(self.flux_a,1e10) # convert to cgs
		self._lambda = np.divide(self.lambda_a,1e10) # convert to meters
		# validate resampled flux
		blanks = 0
		for i in range(self.flux.shape[0]):
			if (self.flux[i] == 0):
				blanks += 1
			elif math.isnan(self.flux[i]):
				self.flux[i] = 0
				blanks += 1
		if (blanks > ((self.flux.shape[0])/5)): # 5th or more values as zero?
			zero_per = blanks / len(self.flux) * 100 # percent of values that are zero
			print('[ error ] : Due to template limitations, {}% of spectrum in this bandpass is zero. Change redshift and/or bandpass.'.format(zero_per))
		else:
			fluxtrans = np.dot(math.fsum(self.flux),self.trans)
			extinction_l = np.dot(math.fsum(self.extinction),math.fsum(self._lambda))
			if (self.mag_sys_opt.lower() == 'vega'):
				vega_panda = pd.read_csv(edl.vega_file,sep='\s+',skiprows=1) # * u.jansky
				flux_vega = np.dot(spectres(self.lambda_a,np.asarray(vega_panda[vega_panda.columns[0]]),np.asarray(vega_panda[vega_panda.columns[1]])),1e10) #.cgs	
				vegatrans = np.dot(math.fsum(self.vega),self.trans)
				self.mag_model = np.add(np.dot(-2.5,np.log10(np.divide(np.dot(fluxtrans,extinction_l)),(np.dot(vegatrans,extinction_l)))),0.026) # assuming that vega is 0.026 mag in all bands
			elif (self.mag_sys_opt.lower() == 'ab'):
				trans_l = np.dot(self.trans,math.fsum(self._lambda))
				om_1 = np.divide(const.c.value,(np.square(self._lambda)))
				om_2 = np.dot(trans_l,math.fsum(self.extinction))
				om_3 = math.fsum(np.divide(np.dot(fluxtrans,extinction_l),np.dot(om_2,om_1)))
				self.mag_model = np.subtract(-48.60,(np.dot(2.5,np.log10(om_3)))) # zeropoint is 48.60
			# compare resampled flux to object's template.
			self.del_mag = self.magnitude - self.mag_model # desired vs template
			self.ss_flux = np.dot(self.flux,(10**(-self.del_mag/2.5))) # rescale
			#except: print('[ error ] : Attribute \'mag_sys_opt\' must be set! Include \'VEGA\' or \'AB\'')
			old_res = self.obj_x[2] - self.obj_x[1]
			sigma = self.delta_lambda / gaussian_sigma_to_fwhm
			self._x = np.arange(start=(-5*sigma),stop=(5*sigma),step=old_res)
			funx_lambda = lambda x,sigma: np.dot((1/(sigma*(math.sqrt(2*math.pi)))),np.exp(np.divide(np.negative(np.square(x)),(2*sigma**2))))
			funx = funx_lambda(self._x,sigma)
			kernel = funx / np.trapz(funx)
			print('[ info ] : Model resolution: {}\n[ info ] : Requested resolution: {}'.format(old_res,self.plot_step))
			if (old_res < self.plot_step):
				self.ss_flux = convolve_fft(self.ss_flux,kernel)
			else:
				self.ss_flux = convolve_fft(self.ss_flux,kernel) # lambda this later for callable response
				print('[ error ] : The object template resolution is lower than the instrument\'s resolution. Reduce the requested resolution (plot_step).')


	def signal(self):
		# counts counts. hehe.
		self.ss_lambda = np.linspace(self.obj_x[0],self.obj_x[-1],self.ss_flux.shape[0]) # interpolate this to match end thing
		ts_flux = (spectres(self.wavelength,self.ss_lambda,self.ss_flux)) # * u.jansky
		power = ((ts_flux * 1e-3) * self.area * self.exposure_time * self.plot_step) # should come out as joules
		self.counts_counts = (power / ((const.h.value * const.c.value) / self.wavelength)) / 1e10
		_sigma = self.seeing / gaussian_sigma_to_fwhm
		lam_func = lambda x: np.dot((1/(_sigma*(math.sqrt(2*math.pi)))),np.exp(np.divide(np.negative(np.square(x)),(2*_sigma**2))))
		self.percent,self.percent_err = np.square(quad(lam_func,(-self.slit_size/2),(self.slit_size/2)))
		print(self.percent) # for signal calculation later
		self.extension = (self.seeing * self.slit_size)
		#mirror_bounce = apply 2
		# final step
		if (self.plot_channel.lower() == 'red'):			
			self.red_total_efficiency = np.dot(np.dot(np.dot(self.red_dichro,self.red_grating),np.dot(self.red_ccd,self.end_extinction)),edl.coating_efficiency)
			self.red_signal = np.dot(np.dot(self.counts_counts,self.percent),self.red_total_efficiency)
		elif (self.plot_channel.lower() == 'blue'):
			self.blue_total_efficiency = np.dot(np.dot(np.dot(self.blue_dichro,self.blue_grating),np.dot(self.blue_ccd,self.end_extinction)),edl.coating_efficiency)
			self.blue_signal = np.dot(np.dot(self.counts_counts,self.percent),self.blue_total_efficiency)
		else:
			self.red_total_efficiency = np.dot(np.dot(np.dot(self.red_dichro,self.red_grating),np.dot(self.red_ccd,self.end_extinction)),edl.coating_efficiency)
			self.red_signal = np.dot(np.dot(self.counts_counts,self.percent),self.red_total_efficiency)
			self.blue_total_efficiency = np.dot(np.dot(np.dot(self.blue_dichro,self.blue_grating),np.dot(self.blue_ccd,self.end_extinction)),edl.coating_efficiency)
			self.blue_signal = np.dot(np.dot(self.counts_counts,self.percent),self.blue_total_efficiency)



	def noise(self):
		# calculate noise
		sigma = self.delta_lambda / gaussian_sigma_to_fwhm
		_x = np.arange(start=(-5*sigma),stop=(5*sigma),step=self.old_res_sky)
		funx_lambda = lambda x,sigma: np.dot((1/(sigma*(math.sqrt(2*math.pi)))),np.exp(np.divide(np.negative(np.square(x)),(2*sigma**2))))
		funx = funx_lambda(_x,sigma)
		kernel = funx / np.trapz(funx)
		#self.sky_flux = convolve(self.sky_y,kernel,'fill') # * u.jansky
		self.sky_flux = spectres(kernel,self.sky_x,self.sky_y) # * u.jansky
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
		self.atmoext()
		self.mirrors()
		self.resample_dichroic()
		self.resample_grating()
		self.resample_ccd()
		self.readnoise()
		self.resample_flux()
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

'''