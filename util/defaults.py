#!/usr/bin/python3
'''
	GMACS ETC: Defaults
		some default values
			if you hadn't picked
				up on that just yet
'''
from numpy import arange

# widget population
string_pagename = 'etc.html'	# filename.html
site_title = 'GMACS ETC'		# like, the webpage 'title' field thing
string_header_one = 'GMACS : Exposure Time Calculator'
string_header_two = 'Munnerlyn Astronomical Instrumentation Lab'
string_title = ['Telescope size','Object type','Grating','Redshift','Exposure time (seconds)','Seeing (arcseconds)',
	'Slit width (arcseconds)','Days from new moon','Wavelength Range (Angstrom)','Binned pixel scale','Calculate',
	'Plot','Galaxy Types']
string_widget_labels = ['Filter']
string_suffixes = ['magnitude','seconds','arcseconds','Angstrom']
string_telescope_sizes = ['Full size (7 mirrors)','First light (4 mirrors)']
string_object_types = ['Stars','Galaxies']
string_object_types_types = ['Star Classifications','Galactic Classifications']
string_star_types = [('O5V','O5V'),('B0V','B0V'),('B5V','B5V'),('A0V','A0V'),('A5V','A5V'),('F0V','F0V'),('F5V','F5V'),
	('G0V','G0V'),('G5V','G5V'),('K0V','K0V'),('K5V','K5V'),('M0V','M0V'),('M5V','M5V')]
string_galaxy_types = [('Starbursts with E(B-V) < 0.10','0'),
	('Starbursts with 0.11 < E(B-V) < 0.21','1'),
	('Starbursts with 0.39 < E(B-V) < 0.50','2'),
	('Starbursts with 0.51 < E(B-V) < 0.60','3'),
	('Starbursts with 0.61 < E(B-V) < 0.70','4'),
	('S0 Galaxies','5'),
	('Sa Galaxies','6'),
	('Sb Galaxies','7'),
	('Sc Galaxies','8'),
	('Galactic Bulges','9'),
	('Elliptical Galaxies','10'),
	('Lyman-break Galaxies','11')]
string_filters = ['u','g','r','i','z','U','B','V','R','I']
string_filters_menu = [('u','0'),('g','1'),('r','2'),('i','3'),('z','4'),('U','5'),('B','6'),('V','7'),('R','8'),('I','9')]
string_magnitude_three = ['Vega','AB']
string_grating_types  = ['Low resolution','High resolution']
string_moon_days = ['0','3','7','10','14']
string_binned_pixel_scale_modes = ['Default','Custom']
string_calculate_types = ['Signal-to-Noise Ratio','Observed Spectrum (without noise)','Observed Spectrum (with noise)','Observed sky background','Dichroic throughput','Grating throughput','CCD QE','Atmospheric extinction']
string_plot_types = ['Red','Blue']
string_binned_pixel_scale_header = 'Binned pixel scale'
string_binned_pixel_scale_labels = ['1x1','2x2']
string_plot_step = 'Plot-Step'
string_axis_labels = [('x','y'),
						('x','y'),
						('x','y'),
						('x','y'),
						('x','y'),
						('x','y'),
						('x','y'),
						('Throughput','Angstroms')] # atmospheric extinction

# widget defaults -- some refer to position in list, some are legit numbers
default_telescope_sizes = string_telescope_sizes[0]
default_object_types = 0
default_star_types = string_star_types[3]
default_galaxy_types = string_galaxy_types[11]
default_filters = string_filters[2]
default_magnitude_two = '25'
default_magnitude_three = string_magnitude_three[1]
default_grating_types = string_grating_types[0]
default_redshift = 0.0
default_exposure_time = 3200
default_exposure_time_start = 0
default_exposure_time_end = 3600
default_seeing = '0.5'
default_slit_width = '0.7'
default_start_wavelength = 2800
default_stop_wavelength = 10360
default_binned_pixel_scale = '0.0' # TODO: actual value
default_calculate_types = string_calculate_types[0]
default_plot_types = [0,1]
default_plot_height = 700
default_plot_width = 1380
default_plot_sizing_mode = 'fixed'
default_toolbar_width = 500
default_plot_step=[0.1,0.01,1]
default_plot_step_step=0.1
default_wavelength = arange(default_start_wavelength,default_stop_wavelength,default_plot_step_step)
default_line_width = 2
default_line_color_red = "red"
default_line_color_blue = "blue"
default_line_color_else = "black"