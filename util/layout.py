#!/usr/bin/python3

'''
	GMACS ETC : Layout
		pretty ummmmmmmmmmm
			pretty self explanatory

'''
import os
import math
import time
import numpy as np 
from collections import defaultdict
from scipy import integrate, interpolate
import pandas as pd
from bokeh.io import reset_output, curdoc
from bokeh.layouts import widgetbox, column, row, layout
from bokeh.models import ColumnDataSource, glyphs
from bokeh.models.widgets import Dropdown, RadioButtonGroup, CheckboxButtonGroup, Slider, TextInput, RangeSlider, Paragraph, Panel, Tabs, Div
from bokeh.models.callbacks import CustomJS
from bokeh.plotting import figure, show
from bokeh.resources import CDN
from bokeh.embed import file_html
from random import choice as tickle # gotta make your own fun
from tqdm import tqdm
# preamble
error_flag = 0

''' default values '''
step_size = 0.1 				# scales interpolation steps
site_title = 'GMACS ETC'		# like, the webpage 'title' field thing
line_width = 2					# might want to allow user to change this
string_header_one = 'GMACS : Exposure Time Calculator'
string_header_two = 'Munnerlyn Astronomical Instrumentation Lab'
# widgit population
string_title = ['Telescope size','Object type','Grating','Redshift','Exposure time (seconds)','Seeing (arcseconds)',
	'Slit width (arcseconds)','Days from new moon','Wavelength Range (Angstrom)','Binned pixel scale','Calculate',
	'Plot','Galaxy Types']
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
string_binned_pixel_scale_manual = 'Enter Any Number:'
string_plot_step = 'Plot-Step'

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
default_exposure_time = 3600
default_exposure_time_start = 0
default_exposure_time_end = 3600
default_seeing = '0.5'
default_slit_width = '0.7'
default_start_wavelength = 3500
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
default_wavelength = np.arange(default_start_wavelength,default_stop_wavelength,default_plot_step_step)

# functional variable initialization
functional_telescope_area = 0
functional_object_type = 0
functional_object_x = []
functional_object_y = []
functional_redshift = default_redshift


'''
widget factory
	TODO: check for float values... I guess?

'''
# radio button groups
widget_telescope_sizes = RadioButtonGroup(labels=string_telescope_sizes, active=0)
widget_object_types = RadioButtonGroup(labels=string_object_types, active=default_object_types)
widget_galaxy_type = Dropdown(label=string_title[12],menu=string_galaxy_types)
widget_mag_type = RadioButtonGroup(labels=string_magnitude_three,active=0)
widget_grating_types = RadioButtonGroup(labels=string_grating_types, active=0)
widget_moon_days = RadioButtonGroup(labels=string_moon_days, active=0) # TODO: need title for this, would be string_title[7]
widget_binned_pixel_scale_mode = RadioButtonGroup(labels=string_binned_pixel_scale_modes, active=0)
# dropdown menus
widget_types_types = Dropdown(label=string_object_types_types[0],menu=string_star_types) # TODO: dynamically set these
widget_filters = Dropdown(label='Filter',menu=string_filters_menu,width=100)
# text input
widget_mag_input = TextInput(value=default_magnitude_two,title=string_suffixes[0].title(),width=75)
widget_redshift = TextInput(value=str(default_redshift), title=string_title[3])
widget_seeing = TextInput(value=default_seeing, title=string_title[5])
widget_slit_width = TextInput(value=default_slit_width, title=string_title[6])
widget_binned_pixel_scale = TextInput(value=default_binned_pixel_scale,title=string_binned_pixel_scale_manual) # TODO: hide on default
# sliders
widget_exposure_time = Slider(start=default_exposure_time_start, end=default_exposure_time_end,
	value=default_exposure_time, step=1, title=string_title[4])
widget_wavelengths = RangeSlider(start=default_start_wavelength, end=default_stop_wavelength,
	value=(default_start_wavelength,default_stop_wavelength), step=1, title=string_title[8])
# other widgets
widget_header = Div(text='<h1>'+string_header_one+'</h1><h3>'+string_header_two+'</h3>',width=500,height=70)
widget_plot_types = CheckboxButtonGroup(labels=string_plot_types, active=default_plot_types) # TODO: make double-off == double-on
widget_message_box = Paragraph(text='')
widget_plot_step = Slider(start=default_plot_step[1], end=default_plot_step[2], value=default_plot_step[0], step=default_plot_step_step, title=string_plot_step)
widget_data_frame = 

'''
plot it
- figures are the... yeah if you're reading this, then you know what that means
- tabN 'Panels' are panel spaces that are flipped between
- the 'tabs' variable builds the tabs bar and points them toward their respective panels

'''
# build figures
p0 = figure(plot_width=default_plot_width, plot_height=default_plot_height,sizing_mode=default_plot_sizing_mode)
p1 = figure(plot_width=default_plot_width, plot_height=default_plot_height,sizing_mode=default_plot_sizing_mode)
p2 = figure(plot_width=default_plot_width, plot_height=default_plot_height,sizing_mode=default_plot_sizing_mode)
p3 = figure(plot_width=default_plot_width, plot_height=default_plot_height,sizing_mode=default_plot_sizing_mode)
p4 = figure(plot_width=default_plot_width, plot_height=default_plot_height,sizing_mode=default_plot_sizing_mode)
p5 = figure(plot_width=default_plot_width, plot_height=default_plot_height,sizing_mode=default_plot_sizing_mode)
p6 = figure(plot_width=default_plot_width, plot_height=default_plot_height,sizing_mode=default_plot_sizing_mode)
p7 = figure(plot_width=default_plot_width, plot_height=default_plot_height,sizing_mode=default_plot_sizing_mode)
# p0.circle([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=20, color="navy", alpha=0.5) # just an example for inits

# assemble graph tabs and panels
try:
	tab0 = Panel(child=p0,title=string_calculate_types[0])
	tab1 = Panel(child=p1,title=string_calculate_types[1])
	tab2 = Panel(child=p2,title=string_calculate_types[2])
	tab3 = Panel(child=p3,title=string_calculate_types[3].title())
	tab4 = Panel(child=p4,title=string_calculate_types[4].title())
	tab5 = Panel(child=p5,title=string_calculate_types[5].title())
	tab6 = Panel(child=p6,title=string_calculate_types[6])
	tab7 = Panel(child=p7,title=string_calculate_types[7].title())
	tabs = Tabs(tabs=[tab0,tab1,tab2,tab3,tab4,tab5,tab6,tab7]) # string_title[10]
except:
	print('[Error :: graph assembly]')
	error_flag += 1


'''
interface layout

'''
# get submenu thing here
try:
	widget_group_one = widgetbox(children=[widget_telescope_sizes,widget_object_types,widget_galaxy_type])
	widget_group_two = layout([[widget_mag_input],[widget_filters,widget_mag_type]])
	widget_group_three = widgetbox(children=[widget_grating_types,widget_redshift,widget_exposure_time,widget_seeing,widget_slit_width,widget_moon_days,widget_wavelengths,widget_binned_pixel_scale_mode,widget_binned_pixel_scale,widget_plot_types,widget_plot_step])
	widgets = column(children=[widget_group_one,widget_group_two,widget_group_three],width=default_toolbar_width)
	inputs = row(children=[widgets,tabs],sizing_mode='scale_height')
except:
	print('[Error :: interface objects]')
	error_flag += 1



''' error handling and closure '''
if error_flag >= 1:
	error_form = ' error'
	if error_flag > 1:
		error_form = error_form+'s'
	else:
		pass
	print('\nWarning! ' + str(error_flag) + error_form + ' occurred during the build.' +
		'\nIdentify the conflicts above and consult the user manual for additional assistance.')
else:
	print('Exporting...')
	l = layout([[widget_header],[inputs]])
	curdoc().add_root(l)
	curdoc().title = site_title
	file_html(l, CDN, "GMACS")
	show(l)
	CustomJS(code='''console.log('GMACS ETC v2.0');''')
	