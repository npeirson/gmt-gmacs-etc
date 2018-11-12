base_data = dict(sb1=dh.galaxy_sb1,sb2=dh.galaxy_sb2,sb3=dh.galaxy_sb3,sb4=dh.galaxy_sb4,
	sb5=dh.galaxy_sb5,sb6=dh.galaxy_sb6,s0=dh.galaxy_s0,sa=dh.galaxy_sa,sb=dh.galaxy_sb,
	sc=dh.galaxy_sc,bulge=dh.galaxy_bulge,ellipticals=dh.galaxy_ellipticals,m5v=dh.m5v,
	lbg_all_flam=dh.galaxy_lbg_all_flam,o5v=dh.o5v,b0v=dh.b0v,b57v=dh.b57v,a0v=dh.a0v,
	a5v=dh.a5v,f0v= h.f0v,f5v=dh.f5v,g0v=dh.g0v,g5v=dh.g5v,k0v=dh.k0v,k5v=dh.k5v,m0v=dh.m0v,
	ff_photon_ux=dh.ff_photon_ux,ff_photon_b=dh.ff_photon_b,ff_photon_v=dh.ff_photon_v,
	ff_photon_r=ff_photon_r,ff_photon_i=ff_photon_i,ff_u=ff_u,ff_g=ff_g,ff_r=ff_r,ff_i=ff_i,
	ff_z=dh.ff_z,skyfile_00=dh.skyfile_00,skyfile_03=dh.skyfile_03,skyfile_07=dh.skyfile_07,
	skyfile_10=dh.skyfile_10,skyfile_14=dh.skyfile_14,dichroic_file=dh.dichroic_file,
	grating_blue=dh.grating_blue,grating_red=dh.grating_red,ccd_blue=dh.ccd_blue,ccd_red=dh.ccd_red,
	atmo_ext_file=dh.atmo_ext_file,mirror_file=dh.mirror_file)

bokeh_inserts = dict(tabs=tabs,widget_wavelengths=widget_wavelengths,widget_exposure_time=widget_exposure_time,
	widget_object_types=widget_object_types,widget_galaxy_type=widget_galaxy_type,widget_types_types=widget_types_types,widget_mag_type=widget_mag_type,
	widget_seeing=widget_seeing,widget_grating_types=widget_grating_types,widget_slit_width=widget_slit_width,widget_telescope_sizes=widget_telescope_sizes,
	widget_filers=widget_filters,widget_plot_types=widget_plot_types,widget_binned_pixel_scale=widget_binned_pixel_scale,widget_withnoise=widget_withnoise,
	widget_moon_days=widget_moon_days,widget_redshift=widget_redshift,widget_mag_input=widget_mag_input)


# this uses a caller singleton and enumerated caller-association groups to update only necessary data.
# the else statement catches errors and changes from wavelength or channel

import time
import math

import values as edl
import defaults as dfs
import datahandler as dh

import numpy as np
#from spectres import spectres # using urs
from scipy import interpolate, integrate

from astropy import constants as astroconst
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_sigma_to_fwhm

from bokeh.plotting import figure
from bokeh.embed import autoload_static
from bokeh.io import output_file,show,save
from bokeh.models.callbacks import CustomJS
from bokeh.layouts import widgetbox, column, row, layout
from bokeh.models import ColumnDataSource, glyphs, ranges
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
widget_withnoise = RadioButtonGroup(name='withnoise',labels=['off','on'],active=0)

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

cds_blue = ColumnDataSource(dict(xb=[],yb=[]))
cds_red = ColumnDataSource(dict(xr=[],yr=[]))
gly_blue = glyphs.Line(x="xb",y="yb",line_width=dfs.default_line_width,line_color='blue')
gly_red = glyphs.Line(x="xr",y="yr",line_width=dfs.default_line_width,line_color='red')
p0.add_glyph(cds_blue,gly_blue)
p0.add_glyph(cds_red,gly_red)

bokeh_inserts = dict(cds_blue=cds_blue,cds_red=cds_red,tabs=tabs,widget_wavelengths=widget_wavelengths,widget_exposure_time=widget_exposure_time,
	widget_object_types=widget_object_types,widget_galaxy_type=widget_galaxy_type,widget_types_types=widget_types_types,widget_mag_type=widget_mag_type,
	widget_seeing=widget_seeing,widget_grating_types=widget_grating_types,widget_slit_width=widget_slit_width,widget_telescope_sizes=widget_telescope_sizes,
	widget_filers=widget_filters,widget_plot_types=widget_plot_types,widget_binned_pixel_scale=widget_binned_pixel_scale,widget_withnoise=widget_withnoise,
	widget_moon_days=widget_moon_days,widget_redshift=widget_redshift,widget_mag_input=widget_mag_input)

initial_values = dict(telescope_mode=0,wavelength=dfs.default_wavelength,exposure_time=3600,object_type='a5v',
	filter_index=3,mag_sys_opt='ab',magnitude=25,redshift=0,seeing=0.5,slit_size=0.5,moon_days=0,grating_opt=0,
	noise=False,bin_option=edl.bin_options_int[edl.bin_options_default_index],channel='both',sss=False)


''' callback management '''


#carrier_bi = ColumnDataSource(initial_values)
cds_wavelength = ColumnDataSource(dict(x=dfs.default_wavelength))
carrier_iv = ColumnDataSource(dict(data=[0,3600,4,3,1,25,0,0.5,0.5,0,0,False,edl.bin_options_int[edl.bin_options_default_index],'both',False]))
singletons_scalar = ColumnDataSource(dict(scalars=[0,0.,0.,0.,0.,0.,0.,0.,0]))
singletons_arrays = ColumnDataSource(dict(extinction=[],power=[],counts=[],counts_noise=[],lambda_A=[],trans=[],_extinction=[],flux=[],_lambda=[],
											flux_y=[],delta_lambda=[],mag_model=[],flux_A=[],selected_filter=[],sky_flux=[],readnoise=[],mirror=[]))
singletons_matrices = ColumnDataSource(dict(snr_blue=[],snr_red=[],signal_blue=[],signal_red=[],noise_blue=[],noise_red=[],error_blue=[],error_red=[],
											dichro_blue=[],dichro_red=[],ccd_blue=[],ccd_red=[],grating_blue=[],grating_red=[],total_eff_blue=[],
											total_eff_red=[],object_x=[],object_y=[]))

'''

carrier_iv = ColumnDataSource(dict(telescope_mode=0,expoure_time=3600,object_type=4,filter_index=3,mag_sys_opt=1,magnitude=25,redshift=0,seeing=0.5,slit_size=0.5,moon_days=0,grating_opt=0,
	noise=False,bin_option=1,channel=0,sss=0))
'''

def convertable_callback(initial_values=carrier_iv,ss=singletons_scalar,sa=singletons_arrays,sm=singletons_matrices,wavelength=cds_wavelength,
	cds_blue=cds_blue,cds_red=cds_red,tabs=tabs,widget_wavelengths=widget_wavelengths,widget_exposure_time=widget_exposure_time,
	widget_object_types=widget_object_types,widget_galaxy_type=widget_galaxy_type,widget_types_types=widget_types_types,widget_mag_type=widget_mag_type,
	widget_seeing=widget_seeing,widget_grating_types=widget_grating_types,widget_slit_width=widget_slit_width,widget_telescope_sizes=widget_telescope_sizes,
	widget_filers=widget_filters,widget_plot_types=widget_plot_types,widget_binned_pixel_scale=widget_binned_pixel_scale,widget_withnoise=widget_withnoise,
	widget_moon_days=widget_moon_days,widget_redshift=widget_redshift,widget_mag_input=widget_mag_input):
	''' tier 0 '''


	def change(caller,ss=ss,sa=sa,sm=sm): # caller=cb_obj.name,tabs=tabs
		caller = cb_obj.value

		if 0 in widget_plot_types.active:
			acp = 1
		if 1 in widget_plot_types.active:
			acp = 2
		else:	
			ss.data['channel'] = 'both'
		if (acp == 3):
			ss.data['channel'] = 'both'
		elif (acp == 2):
			ss.data['channel'] = 'blue'
		elif (acp == 1):
			ss.data['channel'] = 'red'
		else:
			ValueError('{} Channel allocation error!'.format(edl.string_prefix))

		sa.data.change.emit()


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
				plot_y1 = sm.data['snr_blue']
				plot_y2 = sm.data['snr_red']
			elif (tabs.active == 1): # obs spec no noise
				if (widget_withnoise.active == 1):
					plot_y1 = np.add(sm.data['signal_blue'],sm.data['error_blue'])
					plot_y2 = np.add(sm.data['signal_red'],sm.data['error_red'])
				else:
					plot_y1 = sm.data['signal_blue']
					plot_y2 = sm.data['signal_red']
			elif (tabs.active == 2): # obs sky background
				plot_y1 = sm.data['error_blue']
				plot_y2 = sm.data['error_red']
			else:
				pass # pass to other options

		else:
			if (tabs.active == 3): # dichroic throughput
				plot_y1 = sm.data['dichro_blue']
				plot_y2 = sm.data['dichro_red']
			elif (tabs.active == 4): # grating throughput
				plot_y1 = sm.data['grating_blue']
				plot_y2 = sm.data['grating_red']
			elif (tabs.active == 5): # CCD QE
				plot_y1 = sm.data['ccd_blue']
				plot_y2 = sm.data['ccd_red']
			elif (tabs.active == 6): # atmospheric extinction
				plot_y1 = sa.data['extinction']
				plot_y2 = sa.data['extinction']
			else:
				raise ValueError("{} Invalid active tab: {}".format(edl.string_prefix,tabs.active))

			if (cb_obj.name == 'channel'): # change to actual instance's name attribute, not string
				if 0 in cb_obj.active:
					red_active = True
				if 1 in cb_obj.active:
					blue_active = True
				else:
					red_active = True
					blue_active = True
			if red_active and blue_active:
				ss.data['channel'] = 'both'
			elif blue_active and not red_active:
				ss.data['channel'] = 'blue'
			elif red_active and not blue_active:
				ss.data['channel'] = 'red'

		cds_blue.data['xb'] = cds_wavelength.data['x']
		cds_red.data['xr'] = cds_wavelength.data['x']
		cds_blue.data['yb'] = plot_y1
		cds_red.data['yr'] = plot_y2

		ss.data.change.emit()
		cds_blue.data.change.emit()
		cds_red.data.change.emit()


		def recalculate_snr(caller,ss=ss,sm=sm):
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

			if (ss.data['channel'] == 'blue') or (ss.data['channel'] == 'both'):
				sm.data['snr_blue'] = np.divide(sm.data['signal_blue'],np.sqrt(sm.data['signal_blue'] + sm.data['noise_blue'] + np.square(sa.data['readnoise'])))
			if (ss.data['channel'] == 'red') or (ss.data['channel'] == 'both'):
				sm.data['snr_red'] = np.divide(sm.data['signal_red'],np.sqrt(sm.data['signal_red'] + sm.data['noise_red'] + np.square(sa.data['readnoise'])))

			sm.data.change.emit()


		def recalculate_error(caller,ss=ss,sm=sm):
			if caller in edl.signal_keys:
				self.recalculate_signal(caller)
			if caller in edl.noise_keys:
				self.recalculate_noise(caller)
			if caller in edl.readnoise_keys:
				self.recalculate_readnoise(caller)
			if caller in edl.snr_keys:
				self.recalculate_snr(caller)
			else:
				self.recalculate_signal(caller)
				self.recalculate_noise(caller)
				self.recalculate_readnoise(caller)
				self.recalculate_snr(caller)

			if (ss.data['channel'] == 'blue') or (ss.data['channel'] == 'both'):
				sigma_blue = np.sqrt(sm.data['signal_blue'] + sm.data['noise_blue'] + np.square(sa.data['readnoise']))
				sm.data['error_blue'] = np.random.normal(loc=0, scale=sigma_blue,size=len(sm.data['snr_blue']))
			if (ss.data['channel'] == 'red') or (ss.data['channel'] == 'both'):
				sigma_red = np.sqrt(sm.data['signal_red'] + sm.data['noise_red'] + np.square(sa.data['readnoise']))
				sm.data['error_red'] = np.random.normal(loc=0, scale=sigma_red,size=len(sm.data['snr_red']))

			sm.data.change.emit()


		''' tier 1 '''

		def recalculate_signal(caller,ss=ss,sm=sm):
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

			if (ss.data['channel'] == 'blue') or (ss.data['channel'] == 'both'):
				sm.data['signal_blue'] = np.multiply((sa.data['counts'] * ss.data['data'][5]), sm.data['total_eff_blue'])
			if (ss.data['channel'] == 'red') or (ss.data['channel'] == 'both'):
				sm.data['signal_red'] = np.multiply((sa.data['counts'] * ss.data['data'][5]), sm.data['total_eff_red'])

			sm.data.change.emit()


		def recalculate_noise(caller,ss=ss,sm=sm):
			if caller in edl.counts_noise_keys:
				self.recalculate_counts_noise(caller)
			else:
				self.recalculate_counts_noise(caller)
				self.recalculalte_efficiency_noise(caller)

			if (ss.data['channel'] == 'blue') or (ss.data['channel'] == 'both'):
				sm.data['noise_blue'] = np.multiply(sa.data['counts_noise'],sm.data['total_eff_noise_blue'])
			if (ss.data['channel'] == 'red') or (ss.data['channel'] == 'both'):
				sm.data['noise_red'] = np.multiply(sa.data['counts_noise'],sm.data['total_eff_noise_red'])

			sm.data.change.emit()


		''' tier 2 '''

		def recalculate_counts(caller,sa=sa):
			if caller in edl.counts_keys:
				self.recalculate_flux(caller)
				self.change_telescope_mode(caller)
			else:
				self.recalculate_flux(caller)
				self.change_telescope_mode(caller)
				self.change_plot_step(caller)

			sa.data['power'] = sa.data['flux_y'] * ss.data['data'][0] * widget_exposure_time.value * ss.data['data'][7]
			sa.data['counts'] = np.divide(np.divide(sa.data['power'],np.divide((astroconst.h.value * astroconst.c.value),cds_wavelength.data['x'])),1e10)

			sa.data.change.emit()


		def recalculate_counts_noise(caller,sa=sa):
			if (caller == 'moon_days'):
				self.recalculate_sky_flux(caller)
			if (caller == 'telescope_mode'):
				self.change_telescope_mode(caller)
			if caller in edl.extension_keys:
				self.recalculate_seeing(caller)
			else:
				self.recalculate_sky_flux(caller)
				self.recalculate_seeing(caller)
				self.change_telescope_mode(caller)
				self.change_plot_step(caller)

			sa.data['counts_noise'] = np.multiply(np.multiply(sa.data['sky_flux'],ss.data['data'][6]),(ss.data['data'][0]*widget_exposure_time.value*ss.data['data'][7]))

			sa.data.change.emit()


		def recalculate_percent(caller):
			recalculate_seeing(caller)
			# this is a forwarding function


		def recalculate_efficiency(caller,ss=ss,sm=sm):
			if caller in edl.grating_keys:
				self.recalculate_grating(caller)
			else:
				self.recalculate_atmospheric_extinction(caller)
				self.recalculate_grating(caller)
				self.recalculate_dichroic(caller)
				self.recalculate_ccd(caller)
				self.recalculate_atmospheric_extinction(caller)
				self.recalculate_mirror(caller)

			if (ss.data['channel'] == 'blue') or (ss.data['channel'] == 'both'):
				sm.data['total_eff_blue'] = np.multiply(np.multiply(sm.data['dichro_blue'],sm.data['grating_blue']),np.multiply((sm.data['ccd_blue'] * (edl.coating_eff_blue * sa.data['extinction'])),np.square(sa.data['mirror'])))
			if (ss.data['channel'] == 'red') or (ss.data['channel'] == 'both'):
				sm.data['total_eff_red'] = np.multiply(np.multiply(sm.data['dichro_red'],sm.data['grating_red']),np.multiply((sm.data['ccd_red'] * (edl.coating_eff_red * sa.data['extinction'])),np.square(sa.data['mirror'])))

			sm.data.change.emit()


		def recalculalte_efficiency_noise(caller,ss=ss,sm=sm):
			if caller in edl.grating_keys:
				self.recalculate_grating(caller)
			else:
				self.recalculate_dichroic(caller)
				self.recalculate_grating(caller)
				self.recalculate_ccd(caller)
				self.recalculate_mirror(caller)

			if (ss.data['channel'] == 'blue') or (ss.data['channel'] == 'both'):
				sm.data['total_eff_noise_red'] = np.multiply(np.multiply(sm.data['dichro_blue'],sm.data['grating_blue']),(sm.data['ccd_blue'] * np.square(sa.data['mirror']) * edl.coating_eff_blue))
			if (ss.data['channel'] == 'red') or (ss.data['channel'] == 'both'):
				sm.data['total_eff_noise_blue'] = np.multiply(np.multiply(sm.data['dichro_red'],sm.data['grating_red']),(sm.data['ccd_red'] * np.square(sa.data['mirror']) * edl.coating_eff_red))

			sm.data.change.emit()


		def recalculate_flux(caller,ss=ss,sa=sa):
			if caller in edl.grating_opt_keys:
				self.change_grating_opt(caller)
			if caller in edl.filter_keys:
				self.change_filter(caller)
			if (caller == 'object_type'):
				self.change_object_type(caller)
			else: # callers prolly not necessary but good formality
				self.recalculate_atmospheric_extinction(caller)
				self.change_grating_opt(caller)
				self.change_filter(caller)
				self.change_moon_days(caller)
				self.change_object_type(caller)
				self.change_plot_step(caller)

			# heal identicalities
			sa.data['lambda_A'][0] = sa.data['lambda_A'][0] + ss.data['data'][7]
			sa.data['lambda_A'][-1] = sa.data['lambda_A'][-1] - ss.data['data'][7]

			ftrans = interpolate.interp1d(sa.data['selected_filter'][0],sa.data['selected_filter'][1], kind='cubic')
			sa.data['trans'] = ftrans(sa.data['lambda_A'])

			sa.data['_extinction'] = self.spectres(sa.data['lambda_A'],dh.atmo_ext_x,dh.atmo_ext_y)

			sa.data['flux'] = sa.data['flux_A'] * 1e10
			sa.data['_lambda'] = sa.data['lambda_A'] / 1e10

			# valaidate flux
			num_zeros = 0
			for lux in sa.data['flux']:
				if (lux is None) or (lux is 0):
					num_zeros += 1
					lux = 0

			if (num_zeros >= (sa.data['flux'].shape[0]/5)):
				if (num_zeros == sa.data['flux'].shape[0]):
					print('{} Error: No flux in this bandpass!'.format(edl.string_prefix))
					output_flux = np.zeros(cds_wavelength.data['x'].shape[0])
				else:
					percent_zeros = (num_zeros / sa.data['flux'].shape[0]) * 100
					print('{} Note: {}% of this bandpass has zero flux'.format(edl.string_prefix,percent_zeros))

			change_mag_sys_opt(caller)

			del_mag = float(widget_mag_input.value) - sa.data['mag_model']
			output_flux = np.multiply(sm.data['object_y'],10 ** np.negative(del_mag/2.5))
			old_res = sm.data['object_x'][2] - sm.data['object_x'][1]
			if (old_res < ss.data['data'][7]):
				sa.data['flux_y'] = self.spectres(cds_wavelength.data['x'],sm.data['object_x'],(output_flux*1e-03)) # ergs s-1 cm-2 A-1 to J s-1 m-2 A-1
			else:
				sa.data['flux_y'] = self.spectres(cds_wavelength.data['x'],sm.data['object_x'],(output_flux*1e-03))

			sa.data.change.emit()


		''' tier 3 '''

		def change_grating_opt(caller,sa=sa):
			_active = ''
			if widget_grating_types.active in edl.grating_opt_keys[:2]: # low resolution... or is this backards?
				sa.data['delta_lambda'] = edl.dld[0] * widget_slit_width.value / 0.7
				_active = 'low-resolution'
			elif widget_grating_types.active in edl.grating_opt_keys[3:]:
				sa.data['delta_lambda'] = edl.dld[1] * widget_slit_width.value / 0.7
				_active = 'high-resolution'
			else:
				raise ValueError("{} Invalid grating_opt: {}".format(edl.string_prefix,widget_grating_types.active))
			if not self.sss:
				print("{} Grating changed to {}".format(edl.string_prefix,_active))

			sa.data.change.emit()


		def change_mag_sys_opt(caller,sa=sa):
			#print('_lambda: {}\ntrans: {}\n_extinction: {}\nflux: {}'.format(type(sa.data['_lambda']),type(sa.data['trans']),type(sa.data['_extinction']),type(sa.data['flux']))) # for debugging
			if (widget_mag_type.active == 0): # vega
				flux_vega = self.spectres(cds_wavelength.data['x'],dh.vega[0],dh.vega[1]) * 1e10 # fixed... I hope?
				sa.data['mag_model'] = -2.5 * np.log10(np.divide(math.fsum(sa.data['flux'] * sa.data['_extinction'] * sa.data['_lambda'] * sa.data['trans']),math.fsum(flux_vega * sa.data['trans'] * sa.data['_lambda'] * sa.data['_extinction']))) + 0.03
			elif (widget_mag_type.active == 1): # AB
				sq = np.square(sa.data['_lambda'])
				sa.data['mag_model'] = -48.6 - 2.5 * np.log10(math.fsum(sa.data['flux'] * sa.data['trans'] * sa.data['_extinction'] * sa.data['_lambda']) / math.fsum(sa.data['trans'] * sa.data['_lambda'] * sa.data['_extinction'] * (astroconst.c.value/sq)))
			else:
				raise ValueError("{} Invalid magnitude system option (mag_sys_opt): {}".format(edl.string_prefix,widget_mag_type.active))
			if not sss:
				print("{} Magnitude system changed to {}".format(edl.string_prefix,widget_mag_type.active))

			sa.data.change.emit()


		def change_object_type(caller,sa=sa,sm=sm):
			if widget_object_types.active == 0:
				_obj_num = widget_types_types.active
			elif widget_object_types.active == 1:
				_obj_num = (widget_galaxy_type.active + range(len(edl.stellar_files)))
			else:
				raise ValueError('{} Invalid object type numerical representation: {}'.format(edl.string_prefix,widget_object_types.active))

			if _obj_num in edl.stellar_keys:
				index_of = [i for i,name in enumerate(edl.stellar_keys) if _obj_num == name][0]
				object_data = dh.starfiles[index_of]
			elif _obj_num in edl.galactic_keys:
				index_of = [i for i,name in enumerate(edl.galactic_keys) if _obj_num == name][0]
				object_data = dh.galaxyfiles[index_of]
			else:
				raise ValueError("{} Invalid object type: {}".format(edl.string_prefix,_obj_num))

			sm.data['object_x'] = object_data[0] * (float(widget_redshift.value))
			sm.data['object_y'] = object_data[1]
			sa.data['flux_A'] = spectres(sa.data['lambda_A'],sm.data['object_x'],sm.data['object_y'])
			if not sss:
				print("{} Object type changed to {}".format(edl.string_prefix,_obj_num))

			sm.data.change.emit()
			sa.data.change.emit()


		def change_moon_days(caller,sa=sa):
			if widget_moon_days.actives in edl.moon_days_keys:
				sa.data['sky_background'] = dh.skyfiles[(int(np.where(np.asarray(edl.moon_days_keys)==widget_moon_days.actives)[0]))]
			else:
				raise ValueError('{} Invalid number of days since new moon: {}'.format(edl.string_prefix,widget_moon_days.actives))
			self.recalculate_sky_flux(caller)
			if not sss:
				print("{} Days since new moon changed to {}".format(edl.string_prefix,widget_moon_days.actives))

			sa.data.change.emit()


		def change_telescope_mode(caller,ss=ss):
			if widget_telescope_sizes.active == 1:
				ss.data['data'][0] = edl.area[0]
			elif widget_telescope_sizes.active == 0:
				ss.data['data'][0] = edl.area[1]
			else:
				raise ValueError('{} Invalid telescope mode: {}'.format(edl.string_prefix,widget.telescope_sizes.active))
			if not sss:
				print("{} Telescope mode changed to {} (area: {} m^2)".format(edl.string_prefix,widget_telescope_sizes.active,ss.data['data'][0]))

			ss.data.change.emit()


		def change_filter(caller,sa=sa,ss=ss):
			sa.data['selected_filter'] = dh.filterfiles[widget_filters.active]
			filter_min = min(sa.data['selected_filter'][0])
			filter_max = max(sa.data['selected_filter'][0])

			if (filter_min > cds_wavelength.data['x'][0]):
				lambda_min = filter_min
			elif (filter_min == cds_wavelength.data['x'][0]):
				filter_min = sa.data['selected_filter'][int(np.where(sa.data['selected_filter'][0] > cds_wavelength.data['x'][0])[0])]
			else:
				lambda_min = cds_wavelength.data['x'][0]

			if (filter_max < cds_wavelength.data['x'][-1]):
				lambda_max = filter_max
			elif (filter_max == cds_wavelength.data['x'][-1]):
				filter_max = sa.data['selected_filter'][int(np.where(sa.data['selected_filter'][0] < cds_wavelength.data['x'][-1])[-1])]
			else:
				lambda_max = cds_wavelength.data['x'][-1]

			self.change_plot_step(caller) # i dont remember...
			sa.data['lambda_A'] = np.arange(lambda_min,lambda_max,ss.data['data'][7])
			if not sss:
				_active = edl.filter_files[widget_filters.active]
				print("{} Filter changed to {} ({}--{} nm)".format(edl.string_prefix,_active,sa.data['selected_filter'][0][0],sa.data['selected_filter'][0][-1]))

			sa.data.change.emit()


		def recalculate_seeing(caller,ss=ss):
			_sigma = widget_seeing.value / gaussian_sigma_to_fwhm
			funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
			ss.data['data'][1],ss.data['data'][3] = integrate.quad(funx,(-widget_slit_width.value/2),(widget_slit_width.value/2))
			ss.data['data'][2],ss.data['data'][4] = integrate.quad(funx,(-widget_seeing.value/2),(widget_seeing.value/2))
			ss.data['data'][5] = ss.data['data'][1] * ss.data['data'][2] # can use error if you add it later...
			ss.data['data'][6] = widget_seeing.value * widget_slit_width.value
			if not sss:
				print("{} Seeing changed to {}".format(edl.string_prefix,widget_seeing.value))

			ss.data.change.emit()


		def recalculate_sky_flux(self,caller):
			sky_x = sa.data['sky_background'][0] * 1e4
			sky_y = sa.data['sky_background'][1] / 1e4
			old_res = sky_x[2] - sky_x[1]
			_sigma = sa.data['delta_lambda'] / gaussian_sigma_to_fwhm
			_x = np.arange((-5*_sigma),(5*_sigma),old_res)
			funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.multiply(np.square(_sigma),2))))
			degrade = funx(_x)/np.trapz(funx(_x))
			sky_y = convolve_fft(sky_y,degrade)
			sa.data['sky_flux'] = self.spectres(cds_wavelength.data['x'],sky_x,sky_y)

			sa.data.change.emit()


		def recalculate_dichroic(self,caller):
			if (ss.data['channel'] is 'blue') or (ss.data['channel'] is 'both'):
				fblue_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y1, kind='cubic')
				sm.data['dichro_blue'] = fblue_dichro(cds_wavelength.data['x'])
			if (ss.data['channel'] is 'red') or (ss.data['channel'] is 'both'):
				fred_dichro = interpolate.interp1d(dh.dichroic_x,dh.dichroic_y2, kind='cubic')
				sm.data['dichro_red'] = fred_dichro(cds_wavelength.data['x'])

			sm.data.change.emit()


		def recalculate_grating(self,caller):
			if (ss.data['channel'] is 'blue') or (ss.data['channel'] is 'both'):
				sm.data['grating_blue'] = self.spectres(cds_wavelength.data['x'],(dh.grating1[0]*10),dh.grating1[1])
			if (ss.data['channel'] is 'red') or (ss.data['channel'] is 'both'):
				sm.data['grating_red'] = self.spectres(cds_wavelength.data['x'],(dh.grating2[0]*10),dh.grating2[1])

			sm.data.change.emit()


		def recalculate_ccd(self,caller):
			if (ss.data['channel'] is 'blue') or (ss.data['channel'] is 'both'):
				fblue_ccd = interpolate.interp1d((dh.ccd1[0]*10),dh.ccd1[1], kind='cubic')
				sm.data['ccd_blue'] = fblue_ccd(cds_wavelength.data['x'])
			if (ss.data['channel'] is 'red') or (ss.data['channel'] is 'both'):
				fred_ccd = interpolate.interp1d((dh.ccd2[0]*10),dh.ccd2[1], kind='cubic')
				sm.data['ccd_red'] = fred_ccd(cds_wavelength.data['x'])

			sm.data.change.emit()


		def recalculate_mirror(self,caller):
			fmirror = interpolate.interp1d(dh.mirror_file_x,dh.mirror_file_y, kind='cubic')
			sa.data['mirror'] = fmirror(cds_wavelength.data['x'])

			sa.data.change.emit()


		def recalculate_readnoise(self,caller): # probably not implemented in all necessary places yet...
			spectral_resolution = math.ceil((widget_slit_width.value/(0.7/12))/2)*2 # px (ceil()/2)*2 to round up to next even integer
			spatial_resolution = math.ceil((widget_seeing.value/(0.7/12))/2)*2 # probably have to find another method, because pscript is depriciating :(
			extent = widget_seeing.value * widget_slit_width.value
			npix = extent/(0.7/12)**2
			
			rn = edl.rn_default

			if (widget_binned_pixel_scale.active > 0) and (widget_binned_pixel_scale.active < 5):
				sa.data['readnoise'] = math.ceil(rn * spectral_resolution * spatial_resolution / (widget_binned_pixel_scale.active**2))
				if not sss:
					print('{} Pixel binning: ({}x{})'.format(edl.string_prefix,widget_binned_pixel_scale.active,widget_binned_pixel_scale.active))
					print('{} Extent: {} arcsec^2\n{} num pixels/resel: {} px\n{} spectral resolution: {} px\n{} spatial resolution: {} px'.format(edl.string_prefix,extent,edl.string_prefix,int(math.ceil(npix)),edl.string_prefix,spectral_resolution,edl.string_prefix,spatial_resolution))
			else:
				raise ValueError('{} Invalid pixel binning option: {}'.format(edl.string_prefix,widget_binned_pixel_scale.active))

			sa.data.change.emit()


		def recalculate_atmospheric_extinction(self,caller):
			sa.data['extinction'] = self.spectres(cds_wavelength.data['x'],dh.atmo_ext_x,dh.atmo_ext_y)
			
			sa.data.change.emit()


		def change_plot_step(self,caller):
			if (cds_wavelength.data['x'][0] >= dfs.default_limits_wavelength[0]) and (cds_wavelength.data['x'][-1] <= dfs.default_limits_wavelength[1]):
				ss.data['data'][7] = cds_wavelength.data['x'][2] - cds_wavelength.data['x'][1]
			else: # technically shouldn't happen...
				raise ValueError('{} Invalid wavelength extrema ({}--{}). Must be within {}--{}'.format(edl.string_prefix,wavelength[0],wavelength[-1],dfs.default_limits_wavelength[0],dfs.default_limits_wavelength[1]))

			ss.data.change.emit()


		def num_router(self,_obj):
			try:
				if isinstance(_obj,int):
					return _obj
				elif isinstance(_obj,float):
					return _obj
				else:
					raise ValueError("{} Invalid value: {}".format(edl.string_prefix,_obj))
			except:
				raise TypeError("{} Invalid type: {}".format(edl.string_prefix,type(_obj)))



		def type_router(self,_obj,_keys):
			try:
				if isinstance(_obj,str):
					return int(np.where(np.asarray(_keys)==_obj.lower())[0])
				else:
					return self.num_router(_obj)
			except:
				raise TypeError("{} Invalid type: {}".format(edl.string_prefix,type(_obj)))


		# adaptation of SpectRes, from [Adam Carnall](https://github.com/ACCarnall/SpectRes)
		def spectres(self,new_spec_wavs, old_spec_wavs, spec_fluxes):
		    # Arrays of left-hand sides and widths for the old and new bins
		    spec_lhs = np.zeros(old_spec_wavs.shape[0])
		    spec_widths = np.zeros(old_spec_wavs.shape[0])
		    spec_lhs = np.zeros(old_spec_wavs.shape[0])
		    spec_lhs[0] = old_spec_wavs[0]
		    spec_lhs[0] -= (old_spec_wavs[1] - old_spec_wavs[0])/2
		    spec_widths[-1] = (old_spec_wavs[-1] - old_spec_wavs[-2])
		    spec_lhs[1:] = (old_spec_wavs[1:] + old_spec_wavs[:-1])/2
		    spec_widths[:-1] = spec_lhs[1:] - spec_lhs[:-1]

		    filter_lhs = np.zeros(new_spec_wavs.shape[0]+1)
		    filter_widths = np.zeros(new_spec_wavs.shape[0])
		    filter_lhs[0] = new_spec_wavs[0]
		    filter_lhs[0] -= (new_spec_wavs[1] - new_spec_wavs[0])/2
		    filter_widths[-1] = (new_spec_wavs[-1] - new_spec_wavs[-2])
		    filter_lhs[-1] = new_spec_wavs[-1]
		    filter_lhs[-1] += (new_spec_wavs[-1] - new_spec_wavs[-2])/2
		    filter_lhs[1:-1] = (new_spec_wavs[1:] + new_spec_wavs[:-1])/2
		    filter_widths[:-1] = filter_lhs[1:-1] - filter_lhs[:-2]

		    if filter_lhs[0] < spec_lhs[0] or filter_lhs[-1] > spec_lhs[-1]:
		        raise ValueError("[ SpectRes ] : The new wavelengths specified must fall"
		                         "within the range of the old wavelength values.")

		    # Generate output arrays to be populated
		    res_fluxes = np.zeros(spec_fluxes[0].shape + new_spec_wavs.shape)

		    start = 0
		    stop = 0

		    # Calculate new flux and uncertainty values, loop over new bins
		    for j in range(new_spec_wavs.shape[0]):

		        # Find first old bin which is partially covered by the new bin
		        while spec_lhs[start+1] <= filter_lhs[j]:
		            start += 1

		        # Find last old bin which is partially covered by the new bin
		        while spec_lhs[stop+1] < filter_lhs[j+1]:
		            stop += 1

		        # If new bin is fully within one old bin these are the same
		        if stop == start:
		            res_fluxes[j] = spec_fluxes[start]

		        # Otherwise multiply the first and last old bin widths by P_ij
		        else:
		            start_factor = ((spec_lhs[start+1] - filter_lhs[j])
		                            / (spec_lhs[start+1] - spec_lhs[start]))

		            end_factor = ((filter_lhs[j+1] - spec_lhs[stop])
		                          / (spec_lhs[stop+1] - spec_lhs[stop]))

		            spec_widths[start] *= start_factor
		            spec_widths[stop] *= end_factor

		            # Populate res_fluxes spectrum and uncertainty arrays
		            f_widths = spec_widths[start:stop+1]*spec_fluxes[start:stop+1]
		            res_fluxes[j] = np.sum(f_widths, axis=-1)
		            res_fluxes[j] /= np.sum(spec_widths[start:stop+1])

		            # Put back the old bin widths to their initial values for later use
		            spec_widths[start] /= start_factor
		            spec_widths[stop] /= end_factor

		        return res_fluxes
	        # end of simulate.spectres
	print("{} Call from {}".format('me',cb_obj.name))
	change(cb_obj.name)
	# end of class `simulate`


def plot_types_callback(gly_red=gly_red,gly_blue=gly_blue,tabs=tabs):
	# red glyphs on
	if (0 in cb_obj.active):
		gly_red.line_alpha = 0.5
	# red glyphs off
	elif (0 not in cb_obj.active):
		gly_red.line_alpha = 0.0
	# blue glyphs on
	if (1 in cb_obj.active):
		gly_blue.line_alpha = 0.5
	# blue glyphs off
	elif (1 not in cb_obj.active):
		gly_blue.line_alpha = 0.0


''' callback ligatures '''

#widget_plot_types.js_on_change('value',CustomJS.from_py_func(plot_types_callback),args=initial_values)
general_callback = CustomJS.from_py_func(convertable_callback)
widget_seeing.js_on_change('value',general_callback)

widget_withnoise.js_on_change('active',general_callback)

''' final ui panel building '''

widget_group_one = widgetbox(children=[widget_telescope_sizes,widget_object_types,widget_types_types,widget_galaxy_type])
widget_group_two = layout([[widget_mag_input],[widget_filters,widget_mag_type]])
widget_group_three = widgetbox(children=[widget_grating_types,widget_redshift,widget_exposure_time,widget_seeing,widget_slit_width,widget_moon_days_header,widget_moon_days,widget_wavelengths,widget_binned_pixel_scale,widget_plot_types,widget_withnoise]) # removed widget_plot_step and widget_binned_pixel_scale_mode
widgets = column(children=[widget_group_one,widget_group_two,widget_group_three],width=dfs.default_toolbar_width)
inputs = row(children=[widgets,tabs],sizing_mode='scale_height')

l = layout([[widget_header],[inputs]])
show(l)

print('Build completed in {} seconds.'.format(time.time()-time_start))


'''
initial_values = dict(telescope_mode='first',wavelength=dfs.default_wavelength,exposure_time=3600,object_type='a5v',
		filter_index=3,mag_sys_opt='ab',magnitude=25,redshift=0,seeing=0.5,slit_size=0.5,moon_days=0,grating_opt=0,
		noise=False,bin_option=edl.bin_options_int[edl.bin_options_default_index],channel='both',sss=True)
'''

'''
def plot_types_callback(gly_snr_red=gly_snr_red,gly_os_noise_red=gly_os_noise_red,gly_os_nonoise_red=gly_os_nonoise_red,
	gly_sky_red=gly_sky_red,gly_dichroic_red=gly_dichroic_red,gly_grating_red=gly_grating_red,gly_ccd_red=gly_ccd_red,
	gly_snr_blue=gly_snr_blue,gly_os_noise_blue=gly_os_noise_blue,gly_os_nonoise_blue=gly_os_nonoise_blue,
	gly_sky_blue=gly_sky_blue,gly_dichroic_blue=gly_dichroic_blue,gly_grating_blue=gly_grating_blue,gly_ccd_blue=gly_ccd_blue,tabs=tabs):

	red_glyphs = [gly_snr_red,gly_os_noise_red,gly_os_nonoise_red,gly_sky_red,gly_dichroic_red,gly_grating_red,gly_ccd_red]
	blue_glyphs = [gly_snr_blue,gly_os_noise_blue,gly_os_nonoise_blue,gly_sky_blue,gly_dichroic_blue,gly_grating_blue,gly_ccd_blue])
	# red glyphs on
	if (0 in cb_obj.active):
		print(cb_obj.active)
		red_glyphs[tabs.active].line_alpha = 0.5
	# red glyphs off
	elif (0 not in cb_obj.active):
		red_glyphs[tabs.active].line_alpha = 0.0
	# blue glyphs on
	if (1 in cb_obj.active):
		blue_glyphs[tabs.active].line_alpha = 0.5
	# blue glyphs off
	elif (1 not in cb_obj.active):
		blue_glyphs[tabs.active].line_alpha = 0.0
'''