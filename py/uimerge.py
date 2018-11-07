# this uses a caller singleton and enumerated caller-association groups to update only necessary data.
# the else statement catches errors and changes from wavelength or channel

import time
import math

import values as edl
import defaults as dfs
import datahandler as dh

import numpy as np
from spectres import spectres
from scipy import interpolate, integrate

from astropy import constants as const
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_sigma_to_fwhm

from bokeh.models import ColumnDataSource, glyphs
from bokeh.models.widgets import Dropdown, RadioButtonGroup, CheckboxButtonGroup, Slider, TextInput, RangeSlider, Paragraph, Panel, Tabs, Div

''' build user interface '''
print('Building...')
time_start = time.time()

output_file(dfs.string_pagename, title=dfs.site_title, mode='cdn')


''' widgets '''
# radio button groups
widget_telescope_sizes = RadioButtonGroup(labels=dfs.string_telescope_sizes, active=0, name=dfs.string_widget_names[0])
widget_object_types = RadioButtonGroup(labels=dfs.string_object_types, active=dfs.default_object_types, name=dfs.string_widget_names[1])
widget_galaxy_type = Dropdown(label=dfs.string_object_types_types[1],menu=dfs.string_galaxy_types, name=dfs.string_widget_names[2])
widget_galaxy_type.disabled = True # change if galaxies default
widget_mag_type = RadioButtonGroup(labels=dfs.string_magnitude_three,active=0, name=dfs.string_widget_names[3])
widget_grating_types = RadioButtonGroup(labels=dfs.string_grating_types, active=0, name=dfs.string_widget_names[4])
widget_moon_days = RadioButtonGroup(labels=dfs.string_moon_days, active=0, name=dfs.string_widget_names[5]) # TODO: need title for this, would be string_title[7]
widget_binned_pixel_scale_mode = RadioButtonGroup(labels=dfs.string_binned_pixel_scale_modes, active=0)
widget_binned_pixel_scale = RadioButtonGroup(name=dfs.string_binned_pixel_scale_header, labels=dfs.string_binned_pixel_scale_labels,active=dfs.default_binned_pixel_scale_uhh)
# dropdown menus... there's one up there, above, too... could probably move it down here
widget_types_types = Dropdown(label=dfs.string_object_types_types[0],menu=dfs.string_star_types, name=dfs.string_widget_names[6]) # TODO: dynamically set these
widget_filters = Dropdown(label=dfs.string_widget_labels[0],menu=dfs.string_filters_menu,width=100, name=dfs.string_widget_names[7])
# text input
widget_mag_input = TextInput(value=dfs.default_magnitude_two,title=dfs.string_suffixes[0].title(),width=100)
widget_redshift = TextInput(value=str(dfs.default_redshift), title=dfs.string_title[3])
#widget_binned_pixel_scale = TextInput(value=dfs.default_binned_pixel_scale,title=dfs.string_binned_pixel_scale_manual) # TODO: hide on default
# sliders
widget_exposure_time = Slider(start=dfs.default_exposure_time_start, end=dfs.default_exposure_time_end,
	value=dfs.default_exposure_time, step=dfs.default_exposure_time_step, title=dfs.string_title[4], name=dfs.string_widget_names[8])
widget_seeing = Slider(start=dfs.default_seeing_start,end=dfs.default_seeing_end,
	value=dfs.default_seeing, step=dfs.default_seeing_step,title=dfs.string_title[5], name=dfs.string_widget_names[9])
widget_slit_width = Slider(start=dfs.default_slit_width_start,end=dfs.default_slit_width_end,
	value=dfs.default_slit_width,step=dfs.default_slit_width_step,title=dfs.string_title[6], name=dfs.string_widget_names[10])
# range sliders
widget_wavelengths = RangeSlider(start=dfs.default_limits_wavelength[0], end=dfs.default_limits_wavelength[1],
	value=(dfs.default_limits_wavelength), step=dfs.default_wavelength_step,
	title=dfs.string_title[8], name=dfs.string_widget_names[11])

# other widgets
widget_header = Div(text='<h1>'+dfs.string_header_one+'</h1><h3>'+dfs.string_header_two+'</h3>',width=500,height=70)
widget_plot_types = CheckboxButtonGroup(labels=dfs.string_plot_types, active=dfs.default_plot_types, name=dfs.string_widget_names[12]) # TODO: make double-off == double-on
widget_message_box = Paragraph(text='hello')
widget_plot_step = Slider(start=dfs.default_plot_step[1], end=dfs.default_plot_step[2], value=dfs.default_plot_step[0], step=dfs.default_plot_step_step, title=dfs.string_plot_step)
widget_moon_days_header = Paragraph(text=dfs.string_title[7])


''' figures and glyphs '''
# build figures
p0 = figure(plot_width=dfs.default_plot_width, plot_height=dfs.default_plot_height,sizing_mode=dfs.default_plot_sizing_mode,
			x_axis_label=dfs.string_axis_labels[0][0],y_axis_label=dfs.string_axis_labels[0][1])
p1 = figure(plot_width=dfs.default_plot_width, plot_height=dfs.default_plot_height,sizing_mode=dfs.default_plot_sizing_mode,
			x_axis_label=dfs.string_axis_labels[1][0],y_axis_label=dfs.string_axis_labels[1][1])
p2 = figure(plot_width=dfs.default_plot_width, plot_height=dfs.default_plot_height,sizing_mode=dfs.default_plot_sizing_mode,
			x_axis_label=dfs.string_axis_labels[2][0],y_axis_label=dfs.string_axis_labels[2][1])
p3 = figure(plot_width=dfs.default_plot_width, plot_height=dfs.default_plot_height,sizing_mode=dfs.default_plot_sizing_mode,
			x_axis_label=dfs.string_axis_labels[3][0],y_axis_label=dfs.string_axis_labels[3][1])
p4 = figure(plot_width=dfs.default_plot_width, plot_height=dfs.default_plot_height,sizing_mode=dfs.default_plot_sizing_mode,
			x_axis_label=dfs.string_axis_labels[4][0],y_axis_label=dfs.string_axis_labels[4][1])
p5 = figure(plot_width=dfs.default_plot_width, plot_height=dfs.default_plot_height,sizing_mode=dfs.default_plot_sizing_mode,
			x_axis_label=dfs.string_axis_labels[5][0],y_axis_label=dfs.string_axis_labels[5][1])
p6 = figure(plot_width=dfs.default_plot_width, plot_height=dfs.default_plot_height,sizing_mode=dfs.default_plot_sizing_mode,
			x_axis_label=dfs.string_axis_labels[6][0],y_axis_label=dfs.string_axis_labels[6][1])
p7 = figure(plot_width=dfs.default_plot_width, plot_height=dfs.default_plot_height,sizing_mode=dfs.default_plot_sizing_mode,
			x_axis_label=dfs.string_axis_labels[7][0],y_axis_label=dfs.string_axis_labels[7][1],x_range=ranges.Range1d(start=dfs.default_start_wavelength, end=dfs.default_stop_wavelength))
figures = [p0,p1,p2,p3,p4,p5,p6,p7]


# assemble graph tabs
tab0 = Panel(child=p0,title=dfs.string_calculate_types[0])
tab1 = Panel(child=p1,title=dfs.string_calculate_types[1])
tab2 = Panel(child=p2,title=dfs.string_calculate_types[2])
tab3 = Panel(child=p3,title=dfs.string_calculate_types[3].title())
tab4 = Panel(child=p4,title=dfs.string_calculate_types[4].title())
tab5 = Panel(child=p5,title=dfs.string_calculate_types[5].title())
tab6 = Panel(child=p6,title=dfs.string_calculate_types[6])
tab7 = Panel(child=p7,title=dfs.string_calculate_types[7].title())
tabs = Tabs(tabs=[tab0,tab1,tab2,tab3,tab4,tab5,tab6,tab7],name='tabs') # string_title[10]

# snr
cds_snr_red = ColumnDataSource(dict(xr=dfs.default_wavelength,yr=panda_14d[panda_14d.columns[0]]))
cds_snr_blue = ColumnDataSource(dict(xb=dfs.default_wavelength,yb=panda_14d[panda_14d.columns[1]]))
gly_snr_red = glyphs.Line(x="xr",y="yr",line_color='red')
gly_snr_blue = glyphs.Line(x="xb",y="yb",line_color='blue')
p0.add_glyph(cds_snr_red,gly_snr_red)
p0.add_glyph(cds_snr_blue,gly_snr_blue)

# obs spec noise TODO put these on one plot and replace p3 ith manifold
cds_os_noise_red = ColumnDataSource(dict(xr=[],yr=[]))
cds_os_noise_blue = ColumnDataSource(dict(xb=[],yb=[]))
gly_os_noise_red = glyphs.Line(x="xr",y="yr")
gly_os_noise_blue = glyphs.Line(x="xb",y="yb")
p1.add_glyph(cds_os_noise_red,gly_os_noise_red)
p1.add_glyph(cds_os_noise_blue,gly_os_noise_blue)

# obs spec no noise
cds_os_nonoise_red = ColumnDataSource(dict(xr=[],yr=[]))
cds_os_nonoise_blue = ColumnDataSource(dict(xb=[],yb=[]))
gly_os_nonoise_red = glyphs.Line(x="xr",y="yr")
gly_os_nonoise_blue = glyphs.Line(x="xb",y="yb")
p2.add_glyph(cds_os_nonoise_red,gly_os_nonoise_red)
p2.add_glyph(cds_os_nonoise_blue,gly_os_nonoise_blue)

cds_sky_red = ColumnDataSource(dict(xr=[],yr=[]))
cds_sky_blue = ColumnDataSource(dict(xb=[],yb=[]))
gly_sky_red = glyphs.Line(x="xr",y="yr")
gly_sky_blue = glyphs.Line(x="xb",y="yb")
p3.add_glyph(cds_sky_red,gly_sky_red)
p3.add_glyph(cds_sky_blue,gly_sky_blue)

# dichroic throughput
cds_dichroic_red = ColumnDataSource(dict(xr=dichro_x,yr=dichro_y1))
cds_dichroic_blue = ColumnDataSource(dict(xb=dichro_x,yb=dichro_y2))
gly_dichroic_red = glyphs.Line(x="xr",y="yr",line_width=dfs.default_line_width,line_color=dfs.default_line_color_red)
gly_dichroic_blue = glyphs.Line(x="xb",y="yb",line_width=dfs.default_line_width,line_color=dfs.default_line_color_blue)
p4.add_glyph(cds_dichroic_red,gly_dichroic_red)
p4.add_glyph(cds_dichroic_blue,gly_dichroic_blue)

# grating throughput
cds_grating_red = ColumnDataSource(dict(xr=grating_red_1,yr=grating_red_2))
cds_grating_blue = ColumnDataSource(dict(xb=grating_blue_1,yb=grating_blue_2))
gly_grating_red = glyphs.Line(x="xr",y="yr",line_width=dfs.default_line_width,line_color=dfs.default_line_color_red)
gly_grating_blue = glyphs.Line(x="xb",y="yb",line_width=dfs.default_line_width,line_color=dfs.default_line_color_blue)
p5.add_glyph(cds_grating_red,gly_grating_red)
p5.add_glyph(cds_grating_blue,gly_grating_blue)

# ccd efficiency
cds_ccd_red = ColumnDataSource(dict(xr=ccd_efficiency_red_1,yr=ccd_efficiency_red_2))
cds_ccd_blue = ColumnDataSource(dict(xb=ccd_efficiency_blue_1,yb=ccd_efficiency_blue_2))
gly_ccd_red = glyphs.Line(x="xr",y="yr",line_width=dfs.default_line_width,line_color=dfs.default_line_color_red)
gly_ccd_blue = glyphs.Line(x="xb",y="yb",line_width=dfs.default_line_width,line_color=dfs.default_line_color_blue)
p6.add_glyph(cds_ccd_red,gly_ccd_red)
p6.add_glyph(cds_ccd_blue,gly_ccd_blue)

# atmospheric extinction
cds_atmo_ext = ColumnDataSource(dict(x=atmo_ext_x, y=atmo_ext_y))
gly_atmo_ext = glyphs.Line(x="x",y="y",line_width=dfs.default_line_width,line_color=dfs.default_line_color_else)
p7.add_glyph(cds_atmo_ext,gly_atmo_ext)

# other sources for computation pipeline
cds_signal = ColumnDataSource(dict(signal=[]))
cds_noise = ColumnDataSource(dict(noise=[]))


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



class simulate:
	"""
	Class for simulating observatations using GMACS, designed to serve as the back-end to a static bokeh site.

	Parameters
	----------


	"""

	def __init__(self,telescope_mode='first',wavelength=dfs.default_wavelength,exposure_time=3600,object_type='a5v',
		filter_index=3,mag_sys_opt='ab',magnitude=25,redshift=0,seeing=0.5,slit_size=0.5,moon_days=0,grating_opt=0,
		noise=False,bin_option=edl.bin_options_int[edl.bin_options_default_index],channel='both',sss=True,**kwargs):
		
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
		self.__dict__ = dict([arg for arg in kwargs if arg in edl.keys]) # pick up any optionals

		self.change('wavelength') # initializing plots is as easy as this


	''' tier 0 '''

	def change(self): # caller=cb_obj.name,tabs=tabs
		# direct the changed values
		if (tabs.active in [_opt for _opt in range(3)]):
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
				raise ValueError("{} Invalid active tab: {}".format(edl.string_prefix,tabs.active))

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
		if caller in counts_keys:
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
		if caller in edl.atmo_ext_keys:
			self.recalculate_atmospheric_extinction(caller)
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
		if caller in edl.mag_sys_keys():
			self.change_mag_sys_opt(caller)
		if caller in edl.atmo_ext_keys:
			self.change_moon_days(caller)
		else:
			self.change_grating_opt()
			self.change_filter()
			self.change_mag_sys_opt()
			self.change_moon_days()
			self.change_plot_step(caller)

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
			raise ValueError("{} Invalid grating_opt: {}".format(edl.string_prefix,self.grating_opt))
		if not sss:
			print("{} Grating changed to {}".format(edl.string_prefix,_active))


	def change_mag_sys_opt(self,caller):
		if (self.mag_sys_opt == 'vega'):
			flux_vega = spectres(self.wavelength,dh.vega[0],dh.vega[1]) * 1e10 # fixed... I hope?
			self.mag_model = -2.5 * np.log10(np.divide(math.fsum(flux * _extinction * _lambda * trans),math.fsum(flux_vega * trans * _lambda * _extinction))) + 0.03
		elif (self.mag_sys_opt == 'ab'):
			self.mag_model = -48.6 - 2.5 * np.log10(math.fsum(flux * trans * _extinction *_lambda) / math.fsum(trans * _lambda * _extinction * (const.c.value/np.square(_lambda))))
		else:
			raise ValueError("{} Invalid magnitude system option (mag_sys_opt): {}".format(edl.string_prefix,self.mag_sys_opt))
		if not sss:
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
		if not sss:
			print("{} Object type changed to {}".format(edl.string_prefix,self.object_type))


	def change_moon_days(self,caller):
		if self.moon_days in edl.moon_days_keys:
			self.sky_background = dh.skyfiles[(int(np.where(np.asarray(edl.moon_days_keys)==self.moon_days)[0]))]
		else:
			raise ValueError('{} Invalid number of days since new moon: {}'.format(edl.string_prefix,self.moon_days))
		self.recalculate_sky_flux(caller)
		if not sss:
			print("{} Days since new moon changed to {}".format(edl.string_prefix,self.moon_days))


	def change_telescope_mode(self,caller):
		if self.telescope_mode in edl.telescope_mode_keys[:3]:
			self.area = edl.area[0]
		elif self.telescope_mode in edl.telescope_mode_keys[4:]:
			self.area = edl.area[1]
		else:
			raise ValueError('{} Invalid telescope mode: {}'.format(edl.string_prefix,self.telescope_mode))
		if not sss:
			print("{} Telescope mode changed to {} (area: {} m^2)".format(edl.string_prefix,self.telescope_mode,self.area))


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

		self.change_plot_step(caller)
		self.lambda_A = np.arange(lambda_min,lambda_max,self.plot_step)
		if not sss:
			_active = edl.filter_files[self.filter_index]
			print("{} Filter changed to {} ({}--{} nm)".format(edl.string_prefix,_active,selected_filter[0],selected_filter[-1]))


	def recalculate_seeing(self,caller):
		_sigma = self.seeing / gaussian_sigma_to_fwhm
		funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
		self.percent_u,self.percent_err_u = integrate.quad(funx,(-self.slit_size/2),(self.slit_size/2))
		self.percent_l,self.percent_err_l = integrate.quad(funx,(-self.seeing/2),(self.seeing/2))
		self.percent = self.percent_u * self.percent_l # can use error if you add it later...
		self.extension = self.seeing * self.slit_size
		if not sss:
			print("{} Seeing changed to {}".format(edl.string_prefix,self.seeing))


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
			if not sss:
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