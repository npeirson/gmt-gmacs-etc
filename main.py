#!/usr/bin/python3
# bork bork bork

'''
GMACS ETC V2.2 
Exposure Time Calculator for the Giant Magellan Telescope Multi-object Astronomical and Cosmological Spectrograph

Desiged for the Munnerlyn Astronomical Instrumentation Lab at Texas A&M University
Coded by the Kwisatz Haderach in Spring 2018, based on previous versions by Ting Li
There's a user's manual somewhere around here, if you want more info

changelog:	2.0 - python re-design
			2.1 - some plotting functions switched to javascript for easier web integration
		  	2.2 - plotting toolkits changed from matplotlib and scikit to bokeh

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
from random import choice as tickle # gotta make your own fun
from tqdm import tqdm
from core.bork import bork

# hello earth
print('\nGMACS ETC v2.0\n\nGathering fairy-dust...')

'''
file names and paths

'''
time_start = time.time() 		# for elapsed build-time calculations
step_size = 0.1 				# scales interpolation steps
site_title = 'GMACS ETC'		# like, the webpage 'title' field thing
line_width = 2					# might want to allow user to change this
string_header_one = 'GMACS : Exposure Time Calculator'
string_header_two = 'Munnerlyn Astronomical Instrumentation Lab'
galaxy_types_path = 'core/kinney/'
galaxy_types_files = ['SB1','SB2','SB3','SB4','SB5','SB6','S0','Sa','Sb','Sc','bulge','ellipticals','lbg_all_flam']
star_types_path = 'core/pickle/'
star_types_files = ['o5v.dat','b0v.dat','b57v.dat','a0v.dat','a5v.dat','f0v.dat','g0v.dat','g5v.dat','k0v.dat','k5v.dat','m0v.dat','m5v.dat']
filter_path = 'core/filters/'
filter_files = ['photonUX','photonB','photonV','photonR','photonI','u','g','r','i','z']
skyfiles_path = 'core/skybackground/'
skyfiles_files = ['00d','03d','07d','10d','14d']
efficiency_path = 'core/efficiencies/'
efficiency_grating_files = ['grating_red_low_res','grating_blue_low_res']
efficiency_ccd_files = ['e2v_red','e2v_blue']
dichroic_file = 'core/dichroic.txt'
atmospheric_extinction = 'core/atmo_extinction.dat'
lol_units = [' unicorns',' pixies',' bunnies',' puppies',' rainbows',' sparkles',' butterflies']
error_flag = 0


'''
mystical portal of data and parsings!
	
'''
# galaxy types
coalesced_galaxy_types = {}
try:
	for i,galaxy_types_file in tqdm(enumerate(galaxy_types_files),
		desc='Loading object data',ncols=0,unit=tickle(lol_units)):
		_path = os.path.join(galaxy_types_path,"{0}{1}".format(galaxy_types_file,'.txt'))
		_value = pd.read_csv(_path,delimiter=',')
		coalesced_galaxy_types[i] = _value
except:
	print('[Error :: object data]')
	error_flag += 1

# star types
coalesced_star_types = {}
try:
	for i,star_types_file in tqdm(enumerate(star_types_files),
		desc='Loading star data',ncols=0,unit=tickle(lol_units)):
		_path = os.path.join(star_types_path,star_types_file) 
		_value = pd.read_csv(_path,delimiter=',')
		coalesced_star_types[i] = _value
except:
	print('[Error :: star data]')
	error_flag += 1

# filters
coalesced_filters = {}
try:
	for i,filter_file in tqdm(enumerate(filter_files),
		desc='Loading filter data',ncols=0,unit=tickle(lol_units)):
		_path = os.path.join(filter_path,"{0}{1}".format(filter_file,'.txt'))
		_value = pd.read_csv(_path,delimiter=',')
		coalesced_filters[i] = _value
except:
	print('[Error :: filter data]')
	error_flag += 1

# sky files
coalesced_sky_files = {}
try:
	for i,skyfiles_file in tqdm(enumerate(skyfiles_files),
		desc='Loading atmospheric data',ncols=0,unit=tickle(lol_units)):
		_path = os.path.join(skyfiles_path,"{0}{1}".format(skyfiles_file,'.txt'))
		_value = pd.read_csv(_path,delimiter=',')
		coalesced_sky_files[i] = _value
except:
	print('[Error :: atmospheric data]')
	error_flag += 1

# grating efficiency files
coalesced_grating_files = {}
try:
	for i,efficiency_grating_file in tqdm(enumerate(efficiency_grating_files),
		desc='Loading grating efficiencies',ncols=0,unit=tickle(lol_units)):
		_path = os.path.join(efficiency_path,"{0}{1}".format(efficiency_grating_file,'.txt'))
		_value = pd.read_csv(_path,delimiter=',')
		coalesced_grating_files[i] = _value
	# parse before sending to the minions... saves time on callbacks
	grating_red_1 = np.dot(coalesced_grating_files[0]['a'],10)
	grating_red_2 = coalesced_grating_files[0]['b']
	grating_blue_1 = np.dot(coalesced_grating_files[1]['a'],10)
	grating_blue_2 = coalesced_grating_files[1]['b']
except:
	print('[Error :: grating efficiencies]')
	error_flag += 1

# CCD efficiency files
coalesced_ccd_files = {}
try:
	for i,efficiency_ccd_file in tqdm(enumerate(efficiency_ccd_files),
		desc='Loading CCD efficiencies',ncols=0,unit=tickle(lol_units)):
		_path = os.path.join(efficiency_path,"{0}{1}".format(efficiency_ccd_file,'.txt'))
		_value = pd.read_csv(_path,delimiter=',')
		coalesced_ccd_files[i] = _value
	# parse parse parse
	ccd_efficiency_red_1 = np.dot(coalesced_ccd_files[0]['a'],10)
	ccd_efficiency_red_2 = coalesced_ccd_files[0]['b'] # data came in like this
	ccd_efficiency_blue_1 = np.dot(coalesced_ccd_files[1]['a'],10)
	ccd_efficiency_blue_2 = coalesced_ccd_files[1]['d']
except:
	print('[Error :: ccd efficiencies]')
	error_flag += 1

# dichroic efficiency and atmospheric extinction
try:
	# dichro
	coalesced_dichroic = pd.read_csv(dichroic_file,delimiter=',')
	dichro_x = np.dot(coalesced_dichroic['a'],10) # wavelength in A
	dichro_y1 = coalesced_dichroic['b'] # reflectivity, blue channel
	dichro_y2 = coalesced_dichroic['c'] # transmission, red channel
	# atmo ext
	coalesced_atmospheric_extinction = pd.read_csv(atmospheric_extinction,delimiter=',')
	atmo_ext_x = coalesced_atmospheric_extinction['a']
	atmo_ext_y = coalesced_atmospheric_extinction['b']
except:
	print('[Error :: dichroic / atmo. extinction]')
	error_flag += 1


'''
prep some arrays
	(forgive me, ram deities)

'''
# for batch calculations
coalesced_lambda = {}
coalesced_lambda_area = {}
coalesced_flux = {}
coalesced_flux_area = {}
coalesced_ot_interps = {}
coalesced_trans = {}
coalesced_extinction = {}
coalesced_lambda_minimums = []
coalesced_lambda_maximums = []
coalesced_object_x = {}
coalesced_object_y = {}
# functional arrays
fox = []
foy = []


'''
flux machine!
	this calculates the flux, yo

'''
try:
	# plot span and step
	for i,coalesced_filter in tqdm(enumerate(coalesced_filters),
		desc='Coalescing filters',ncols=0,unit=tickle(lol_units)):
		coalesced_lambda_minimums = np.append(coalesced_lambda_minimums, np.min(coalesced_filters[i].a))
		coalesced_lambda_maximums = np.append(coalesced_lambda_maximums, np.max(coalesced_filters[i].a))
		coalesced_lambda_area[i] = np.arange(coalesced_lambda_minimums[i],coalesced_lambda_maximums[i],step_size)
	# divide (n,x,y) into {(n,x),(n,y)}
	for i,coalesced_galaxy_type in tqdm(enumerate(coalesced_galaxy_types),
		desc='Coalescing objects',ncols=0,unit=tickle(lol_units)):
		coalesced_object_x[i] = coalesced_galaxy_types[i].a
		coalesced_object_y[i] = coalesced_galaxy_types[i].b
	# interpolate curves based on the first two for-loops
	for i in tqdm(range(0,len(coalesced_galaxy_types)),
		desc='Calculating flux',ncols=0,unit=tickle(lol_units)):
		coalesced_ot_interps[i] = interpolate.interp1d(x=coalesced_galaxy_types[i].a,y=coalesced_galaxy_types[i].b,kind='linear',fill_value='extrapolate')	
		for j in range(0,len(coalesced_lambda_area)):
			coalesced_flux_area[i] = coalesced_ot_interps[i](coalesced_lambda_area[j])
except:
	print('[Error :: data coalescence]')
	error_flag += 1

for i,coalesced_filter in tqdm(enumerate(coalesced_filters),
	desc='Integrating over lambda',ncols=0,unit=tickle(lol_units)):
	coalesced_trans[i] = interpolate.PchipInterpolator(coalesced_filters[i].a,coalesced_filters[i].b)(coalesced_lambda_area[i])

# linear interpolation
extinction = interpolate.interp1d(x=coalesced_atmospheric_extinction.a,y=coalesced_atmospheric_extinction.b,kind='linear')

for i in tqdm(range(0,len(coalesced_atmospheric_extinction)),
	desc='Calculating atmo extinction',ncols=0,unit=tickle(lol_units)):
	for j in range(0,len(coalesced_lambda_area)):
		coalesced_extinction[i] = extinction(coalesced_lambda_area[j])

# unit conversions
for i in tqdm(range(0,len(coalesced_flux_area)),
	desc='Calculating flux',ncols=0,unit=tickle(lol_units)):
	coalesced_flux[i]= np.dot(coalesced_flux_area[i],(1**10)) # converts to erg s^-1 cm^-2 m^-1
for i in tqdm(range(0,len(coalesced_lambda_area)),
	desc='Converting range',ncols=0,unit=tickle(lol_units)):
	coalesced_lambda[i] = np.divide(coalesced_lambda_area[i],(1**10)) # converts to meters

# correct for NaNs
NaN_counter = np.zeros(len(coalesced_flux))
NaN_error = np.zeros(len(coalesced_flux))
for i in tqdm(range(0,len(coalesced_flux)),
	desc='Cleaning bad sectors',ncols=0,unit=tickle(lol_units)):
	for j in range(0,len(coalesced_flux[i])):
		if (math.isnan(coalesced_flux[i][j]) or coalesced_flux[i][j] == 0.0):
			coalesced_flux[i][j] = 0.0  # due to redshift, the flux exceeds some wavelength ranges, so they go to zero
			NaN_counter[i] += 1
		# could do 'elif flux == 0' but i think that isn't a possible interaction
	if (NaN_counter[i] > np.divide(len(coalesced_flux[i]),5)): # no more than 1/5th NaN content
			NaN_error[i] = 1
				
# from here on out, stuff needs to be handled on-demand
'''
except:
	print('[Error :: flux machine]')
	error_flag += 1
'''

'''
call forth the unicorns!
	this is where the interface stuff begins

'''
moon_day = '00'
delta_lambda_low_default =  3.73
delta_lambda_high_default = 1.40

sigma_div = 2.35482
coating_efficiency = 0.8 * (0.98**14)
planck = 6.626 * 10**(-34) ## TODO: arith
celerity = 29979245

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
	CustomJS(code='''console.log('GMACS ETC v2.0');''')
	

'''
area thing
	i dunno why i think i need to do this seperately like this
	but i do, so i did.

'''
def area_thing():
	if (widget_telescope_sizes.active == 0): # big-kids' telescope
		functional_telescope_area = 368
	else: # lil scope
		functional_telescope_area = 222
	return functional_telescope_area

# same thing, except this is for the Vega/AB mode flux modifiers
flux_vega = {}
mag_model = {}
widget_mag_type_active = widget_mag_type.active
def get_the_flux_out(coalesced_flux,coalesced_trans,widget_mag_type_active,mag_model):
	if (widget_mag_type_active == 0): # Vega mode
		_path = os.path.join('core/','alpha_lyr_stis_005.ascii')
		vega = pd.read_csv(_path,delimiter=',')
		for i,lam in enumerate(coalesced_lambda):
			flux_vega[i] = np.dot(interpolate.PchipInterpolator(vega['a'],vega['b'])(lam),1e10) # assuming that vega is 0.026 mag in all bands
			'''
			mag_model[i] = np.dot(-2.5,math.log((np.divide(np.dot(coalesced_flux[i],coalesced_trans[i])
				.dot(flux_vega[i]).dot(coalesced_trans[i])
				.dot(coalesced_extinction,coalesced_lambda[i])))),10) + 0.03
			'''
			# this is the perp!
			mag_model[i] = np.dot(-2.5, math.log(
				np.divide(
					np.dot(coalesced_flux[i],coalesced_trans[i]).dot(coalesced_extinction).dot(lam).dot(flux_vega[i]),
					10))) + 0.03
	
	elif (widget_mag_type_active == 1): # AB mode
		for i,local in enumerate(coalesced_lambda):
			mag_model[i] = -(48.60)-(np.dot(-2.5, math.log(np.dot(coalesced_flux[i],coalesced_trans[i]).dot(coalesced_extinction).dot(local).dot(np.divide(celerity,np.dot(coalesced_lambda_area[i]))),10)))
			'''
			mag_model[i] = -(48.60)-(2.5)*math.log(np.divide(np.dot(coalesced_flux[i],coalesced_trans[i])
				.dot(coalesced_extinction).dot(local).dot(coalesced_trans[i]).dot(local)
				.dot(coalesced_extinction).dot(np.divide(session.celerity,(np.dot(coalesced_lambda_area[i],coalesced_lambda_area[i]))))),10) # zeropoint is 48.60
			'''	
	else:
		print('[Error] :: NaN correction')
	for i,model in enumerate(mag_model):
		functional_mag_value[i] = widget_mag_input.value - model
	return funtional_mag_value[widget_mag_type_active]


'''
process request
	(explicitly partitioned for readability)
	starts a bork session to host the interactivity
	pretty neat, huh? this is how the gangstas code
	on the west-coast---represent

'''
session = bork(mag_sys_opt=widget_mag_type.active, mag=widget_mag_input.value,
	grating_type=widget_grating_types.active, redshift=widget_redshift.value,
	object_type=widget_object_types.active, filter_opt=widget_filters.value,
	plot_step=widget_plot_step.value,slit_size=widget_slit_width.value,
	start_lambda=widget_wavelengths.value[0],end_lambda=widget_wavelengths.value[1],
	seeing=widget_seeing.value,area=area_thing(),exp_time=widget_exposure_time.value,
	active_tab=tabs.active,wavelength=default_wavelength,colors=widget_plot_types.active,
	sky_files=coalesced_sky_files[widget_moon_days.active])


'''
widget handling
	kinda like oompa-loompas

'''
# function first
def handler(attr,old,new):
	# TODO: add data clears before new data posts
		# signal-to-noise ratio, homie
	session.get_wavelength()
	if (session.active_tab == 0):
		get_the_flux_out(coalesced_flux=coalesced_flux,coalesced_trans=coalesced_trans,widget_mag_type_active=widget_mag_type.active,mag_model=mag_model)
		session.snr()
		if (session.colors == 0 or session.colors == 2):
			snr_red_line = Line(session.wavelength(),session.snr_red(),line_width=line_width)
			p0.add_glyph(session.snr_red,snr_red_line)
		if (session.colors == 1 or session.colors == 2):
			snr_blue_line = Line(session.get_wavelength(),session.get_snr_blue(),line_width=line_width)
			p0.add_glyph(session.snr_blue,snr_blue_line)
	# observed spectrum w/o noise... underscores are just to prevent variable-name mixups
	elif (session.active_tab == 1):
		session.signal()
		if (session.colors == 0 or session.colors == 2):
			_signal_red = session.signal_red
			_signal_red = _signal_red[:,np.newaxis]
			_signal_red_line = Line(session.wavelength,_signal_red,line_width=line_width)
			p1.add_glyph(_signal_red,_signal_red_line)
		if (session.colors == 1 or session.colors == 2):
			_signal_blue = session.signal_blue + session.snr_based_error_blue
			_signal_blue = _signal_blue[:,np.newaxis]
			_signal_blue_line = Line(session.wavelength,_signal_blue,line_width=line_width)
			p1.add_glyph(_signal_blue,_signal_blue_line)
	# observed spectrum with noise
	elif (session.active_tab == 2):
		session.get_snr_based_error()
		if (session.colors == 0 or session.colors == 2):
			oswn_red = session.signal_red + session.snr_based_error_red
			oswn_red_line = Line(session.wavelength,oswn_red,line_width=line_width)
			p2.add_glyph(oswn_red,oswn_red_line)
		if (session.colors == 1 or session.colors == 2):
			oswn_blue = session.signal_blue + session.snr_based_error_blue
			oswn_blue_line = Line(session.wavelength,oswn_blue,line_width=line_width)
			p2.add_glyph(oswn_blue,oswn_blue_line)
	# observed sky background
	elif (session.active_tab == 3):
		session.get_noise()
		if (session.colors() == 0 or session.colors() == 2):
			sky_background_red = session.noise_red
			sky_background_red_line = Line(session.wavelength,noise_red,line_width=line_width)
			p3.add_glyph(sky_background_red,sky_background_red_line)
		if (session.colors() == 1 or session.colors() == 2):
			sky_background_blue = session.noise_blue
			sky_background_blue_line = Line(session.wavelength,blue_noise,line_width=line_width)
			p3.add_glyph(sky_background_blue,sky_background_blue_line)
	# dichroic throughput
	elif (session.active_tab == 4):
		session.get_dichro()
		if (session.colors() == 0 or session.colors() == 2):
			dichro_red = interpolate.PchipInterpolator(dichro_x,dichro_y1)(session.wavelength)
			dichro_red_line = Line(session.wavelength,dichro_red,line_width=line_width)
			p4.add_glyph(dichro_red,dichro_red_line)
		if (session.colors() == 1 or session.colors() == 2):
			dichro_blue = interpolate.PchipInterpolator(dichro_x,dichro_y2)(session.wavelength)
			dichro_blue_line = Line(session.wavelength,dichro_blue,line_width=line_width)
			p4.add_glyph(dichro_blue,dichro_blue_line)
	# grating throughput
	elif (session.active_tab == 5):
		if (session.colors() == 0 or session.colors() == 2):
			grating_red = interpolate.PchipInterpolator(grating_red_1,grating_red_2)(session.wavelength)
			grating_red_line = Line(session.wavelength,grating_red,line_width=line_width)
			p5.add_glyph(grating_red,grating_red_line)
		if (session.colors() == 1 or session.colors() == 2):
			grating_blue = interpolate.PchipInterpolator(grating_blue_1,grating_blue_2)(session.wavelength)
			grating_blue_line = Line(session.wavelength,grating_blue,line_width=line_width)
			p5.add_glyph(grating_blue,grating_blue_line)			
	# CCD QE
	elif (session.active_tab == 6):
		if (session.colors() == 0 or session.colors() == 2):
			ccd_red = interpolate.PchipInterpolator(ccd_efficiency_red_1,ccd_efficiency_red_2)(session.wavelength)
			ccd_red_line = Line(session.wavelength,ccd_red,line_width=line_width)
			p6.add_glyph(ccd_red,ccd_red_line)
		if (session.colors() == 1 or session.colors() == 2):
			ccd_blue = interpolate.PchipInterpolator(ccd_efficiency_blue_1,ccd_efficiency_blue_2)(session.wavelength)
			grating_blue_line = Line(session.wavelength,ccd_blue,line_width=line_width)
			p6.add_glyph(ccd_blue,ccd_blue_line)			
	# atmospheric extinction
	elif (session.active_tab == 7):
		_x = coalesced_extinction[session.object_type]
		p7.line(_x,_x,line_width=line_width)
	else:
		print('[Error] :: main interactivity handler, type 1')


def number_handler(attr,old,new):
	_number_inputs = [widget_mag_input.value,widget_exposure_time.value,
	widget_seeing.value,widget_slit_width.value,widget_binned_pixel_scale.value]
	for thing in _number_inputs:
		if (isinstance(thing, float) != True and isinstance(thing, int) != True): # could probably improve this logic
			widget_message_box.text = 'Error: all manual values must be numeric.'


# linkages
widget_telescope_sizes.on_change("active", handler)
widget_object_types.on_change("active", handler)
widget_mag_type.on_change("active", handler)
widget_mag_input.on_change("value",number_handler)
widget_galaxy_type.on_change("value",handler)
widget_grating_types.on_change("active",handler)
widget_redshift.on_change("value",number_handler)
widget_exposure_time.on_change("value",number_handler)
widget_seeing.on_change("value",number_handler)
widget_slit_width.on_change("value",number_handler)
widget_binned_pixel_scale.on_change("value",number_handler)
widget_moon_days.on_change("active",handler)
tabs.on_change("active",handler)

print('\n[Success!] :: Build completed in ' + str(time.time() - time_start) + ' seconds.')
