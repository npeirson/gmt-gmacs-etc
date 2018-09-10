#!/usr/bin/python3
# bork bork bork

'''
GMACS ETC V2.2 

09/10/18 -- this file is being reworked and will not function in current form

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
