import time
import math
import json
import numpy as np
import values as edl
import datahandler as dh
import defaults as dfs
from spectres import spectres
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
from astropy import constants as astroconst
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_sigma_to_fwhm

from bokeh.layouts import widgetbox, column, row, layout
from bokeh.models import ColumnDataSource, glyphs, ranges
from bokeh.models.widgets import Dropdown, RadioButtonGroup, CheckboxButtonGroup, Slider, TextInput, RangeSlider, Paragraph, Panel, Tabs, Div
from bokeh.plotting import figure
from bokeh.embed import autoload_static
from bokeh.io import output_file,show,save
from bokeh.models.callbacks import CustomJS

print('GMACS Exposure Time Calculator\n\tBuilding...')
time_start = time.time()

# radio button groups
widget_telescope_sizes = RadioButtonGroup(labels=dfs.string_telescope_sizes, active=0, name=dfs.string_widget_names[0])
widget_object_types = RadioButtonGroup(labels=dfs.string_object_types, active=dfs.default_object_types, name=dfs.string_widget_names[1])
widget_galaxy_type = Dropdown(label=dfs.string_object_types_types[1],menu=dfs.string_galaxy_types, name=dfs.string_widget_names[2])
widget_galaxy_type.disabled = True # change if galaxies default
widget_mag_type = RadioButtonGroup(labels=dfs.string_magnitude_three,active=0, name=dfs.string_widget_names[3])
widget_grating_types = RadioButtonGroup(labels=dfs.string_grating_types, active=0, name=dfs.string_widget_names[4])
widget_moon_days = RadioButtonGroup(labels=dfs.string_moon_days, active=0, name=dfs.string_widget_names[5]) # TODO: need title for this, would be string_title[7]
#widget_binned_pixel_scale_mode = RadioButtonGroup(labels=dfs.string_binned_pixel_scale_modes, active=0)
widget_binned_pixel_scale = RadioButtonGroup(name=dfs.string_binned_pixel_scale_header, labels=dfs.string_binned_pixel_scale_labels,active=1)
# dropdown menus... there's one up there, above, too... could probably move it down here
widget_types_types = Dropdown(label=dfs.string_object_types_types[0],menu=dfs.string_star_types, name=dfs.string_widget_names[6]) # TODO: dynamically set these
widget_filters = Dropdown(label=dfs.string_widget_labels[0],menu=dfs.string_filters_menu,width=100, name=dfs.string_widget_names[7])
# text input
widget_mag_input = TextInput(value=dfs.default_magnitude_two,title=dfs.string_suffixes[0].title(),width=100)
widget_redshift = TextInput(value=str(dfs.default_redshift), title=dfs.string_title[3])
widget_withnoise = RadioButtonGroup(name='withnoise',labels=['OFF','ON'],active=0)

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

#wavelength = ranges.Range1d(start=widget_wavelengths.value[0],end=widget_wavelengths.value[1],name='wavelength')

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
'''
p7 = figure(plot_width=dfs.default_plot_width, plot_height=dfs.default_plot_height,sizing_mode=dfs.default_plot_sizing_mode,
			x_axis_label=dfs.string_axis_labels[7][0],y_axis_label=dfs.string_axis_labels[7][1],x_range=ranges.Range1d(start=dfs.default_start_wavelength, end=dfs.default_stop_wavelength))
'''


figures = [p0,p1,p2,p3,p4,p5,p6]

# assemble graph tabs
tab0 = Panel(child=p0,title=dfs.string_calculate_types[0])
tab1 = Panel(child=p1,title=dfs.string_calculate_types[1])
tab2 = Panel(child=p2,title=dfs.string_calculate_types[3])
tab3 = Panel(child=p3,title=dfs.string_calculate_types[4].title())
tab4 = Panel(child=p4,title=dfs.string_calculate_types[5].title())
tab5 = Panel(child=p5,title=dfs.string_calculate_types[6].title())
tab6 = Panel(child=p6,title=dfs.string_calculate_types[7])
tabs = Tabs(tabs=[tab0,tab1,tab2,tab3,tab4,tab5,tab6],name='tabs') # string_title[10]


''' ColumnDataSources '''
wave_range = np.arange(widget_wavelengths.value[0],widget_wavelengths.value[1],50)
cds_wavelength = ColumnDataSource(dict(x=wave_range))
cds_blue = ColumnDataSource(dict(xb=[],yb=[]))
cds_red = ColumnDataSource(dict(xr=[],yr=[]))
singletons_scalar = ColumnDataSource(dict(scalars=[0,0.,0.,0.,0.,0.,0.,0.,0]))
singletons_arrays = ColumnDataSource(dict(extinction=[],power=[],counts=[],counts_noise=[],lambda_A=[],trans=[],_extinction=[],flux=[],_lambda=[],
											flux_y=[],delta_lambda=[],mag_model=[],flux_A=[],selected_filter=[],sky_flux=[],readnoise=[],mirror=[]))
singletons_matrices = ColumnDataSource(dict(snr_blue=[],snr_red=[],signal_blue=[],signal_red=[],noise_blue=[],noise_red=[],error_blue=[],error_red=[],
											dichro_blue=[],dichro_red=[],ccd_blue=[],ccd_red=[],grating_blue=[],grating_red=[],total_eff_blue=[],
											total_eff_red=[],object_x=[],object_y=[]))

# function arbiters
cds_spectres = ColumnDataSource(dict(n=[0,0,0,0,0],o=[0,0,0,0,0,],s=[0,0,0,0,0],f=[0,0,0,0,0]))
cds_trapz = ColumnDataSource(dict(_in=[],_out=[]))


def alt_trapz(cds_trapz=cds_trapz):
	x = cds_trapz.data['_in']
	d = x[2] - x[1]
	slice1 = [slice(None)]
	slice2 = [slice(None)]
	slice1[0] = slice(1, None)
	slice2[0] = slice(None, -1)
	cds_trapz.data['_out'] = math.fsum(d * (x[tuple(slice1)] + x[tuple(slice2)]) / 2.0)
	cds_trapz.change.emit()


# adaptation of SpectRes, from [Adam Carnall](https://github.com/ACCarnall/SpectRes)
def spectres(cds_spectres=cds_spectres):
	new_spec_wavs = cds_spectres.data['n']
	old_spec_wavs = cds_spectres.data['o']
	spec_fluxes = cds_spectres.data['s']

	spec_lhs = [i*0 for i in range(len(old_spec_wavs[0]))]
	spec_widths = [i*0 for i in range(len(old_spec_wavs[0]))]
	spec_lhs[0] = old_spec_wavs[0]
	spec_lhs[0] -= (old_spec_wavs[1] - old_spec_wavs[0])/2
	spec_widths[-1] = (old_spec_wavs[-1] - old_spec_wavs[-2])
	spec_lhs[1:] = (old_spec_wavs[1:] + old_spec_wavs[:-1])/2
	spec_widths[:-1] = spec_lhs[1:] - spec_lhs[:-1]
	filter_lhs = [i*0 for i in range(len(new_spec_wavs[0]+1))]
	filter_widths = [i*0 for i in range(len(new_spec_wavs[0]))]
	filter_lhs[0] = new_spec_wavs[0]
	filter_lhs[0] -= (new_spec_wavs[1] - new_spec_wavs[0])/2
	filter_widths[-1] = (new_spec_wavs[-1] - new_spec_wavs[-2])
	filter_lhs[-1] = new_spec_wavs[-1]
	filter_lhs[-1] += (new_spec_wavs[-1] - new_spec_wavs[-2])/2
	filter_lhs[1:-1] = (new_spec_wavs[1:] + new_spec_wavs[:-1])/2
	filter_widths[:-1] = filter_lhs[1:-1] - filter_lhs[:-2]
	if filter_lhs[0] < spec_lhs[0] or filter_lhs[-1] > spec_lhs[-1]:
		raise ValueError("[ SpectRes ] : The new wavelengths specified must fall within the range of the old wavelength values.")
	res_fluxes = [i*0 for i in range(len(len(spec_fluxes[0]) + len(new_spec_wavs.shape)))]
	start = 0
	stop = 0
	for j in range(len(new_spec_wavs[0])):
		while spec_lhs[start+1] <= filter_lhs[j]:
			start += 1
		while spec_lhs[stop+1] < filter_lhs[j+1]:
			stop += 1
		if stop == start:
			res_fluxes[j] = spec_fluxes[start]
		else:
			start_factor = ((spec_lhs[start+1] - filter_lhs[j]) / (spec_lhs[start+1] - spec_lhs[start]))
			end_factor = ((filter_lhs[j+1] - spec_lhs[stop]) / (spec_lhs[stop+1] - spec_lhs[stop]))
			spec_widths[start] *= start_factor
			spec_widths[stop] *= end_factor
			f_widths = spec_widths[start:stop+1]*spec_fluxes[start:stop+1]
			res_fluxes[j] = math.fsum(f_widths[-1])
			res_fluxes[j] /= math.fsum(spec_widths[start:stop+1])
			spec_widths[start] /= start_factor
			spec_widths[stop] /= end_factor

		cds_spectres.data['f'] = res_fluxes
		cds_spectres.change.emit()
        

alt_spectres = CustomJS.from_py_func(spectres)
alt_trapz = CustomJS.from_py_func(alt_trapz)


def general_callback(cds_spectres=cds_spectres,spectres=alt_spectres,trapz=alt_trapz,sb1=dh.galaxy_sb1,sb2=dh.galaxy_sb2,sb3=dh.galaxy_sb3,sb4=dh.galaxy_sb4,
	sb5=dh.galaxy_sb5,sb6=dh.galaxy_sb6,s0=dh.galaxy_s0,sa=dh.galaxy_sa,sb=dh.galaxy_sb,
	sc=dh.galaxy_sc,bulge=dh.galaxy_bulge,ellipticals=dh.galaxy_ellipticals,m5v=dh.m5v,
	lbg_all_flam=dh.galaxy_lbg_all_flam,o5v=dh.o5v,b0v=dh.b0v,b57v=dh.b57v,a0v=dh.a0v,
	a5v=dh.a5v,f0v=dh.f0v,f5v=dh.f5v,g0v=dh.g0v,g5v=dh.g5v,k0v=dh.k0v,k5v=dh.k5v,m0v=dh.m0v,
	ff_photon_ux=dh.ff_photon_ux,ff_photon_b=dh.ff_photon_b,ff_photon_v=dh.ff_photon_v,
	ff_photon_r=dh.ff_photon_r,ff_photon_i=dh.ff_photon_i,ff_u=dh.ff_u,ff_g=dh.ff_g,ff_r=dh.ff_r,ff_i=dh.ff_i,
	ff_z=dh.ff_z,skyfile_00=dh.skyfile_00,skyfile_03=dh.skyfile_03,skyfile_07=dh.skyfile_07,
	skyfile_10=dh.skyfile_10,skyfile_14=dh.skyfile_14,dichroic_file=dh.dichroic_file,
	grating_blue=dh.grating_blue,grating_red=dh.grating_red,ccd_blue=dh.ccd_blue,ccd_red=dh.ccd_red,
	atmo_ext_file=dh.atmo_ext_file,mirror_file=dh.mirror_file,tabs=tabs,widget_wavelengths=widget_wavelengths,
	widget_exposure_time=widget_exposure_time,widget_object_types=widget_object_types,
	widget_galaxy_type=widget_galaxy_type,widget_types_types=widget_types_types,widget_mag_type=widget_mag_type,
	widget_seeing=widget_seeing,widget_grating_types=widget_grating_types,widget_slit_width=widget_slit_width,
	widget_telescope_sizes=widget_telescope_sizes,widget_filters=widget_filters,widget_plot_types=widget_plot_types,
	widget_binned_pixel_scale=widget_binned_pixel_scale,widget_withnoise=widget_withnoise,
	widget_moon_days=widget_moon_days,widget_redshift=widget_redshift,widget_mag_input=widget_mag_input,
	cds_blue=cds_blue,cds_red=cds_red,cds_wavelength=cds_wavelength):

	string_prefix = '[ etc ] :'
	coating_eff_red = 0.62
	coating_eff_blue = 0.60

	wavelength = range(widget_wavelengths.value[0],widget_wavelengths.value[1],50)
	'''
	wavelength = cds_wavelength.data['x']
	if (wavelength[0] != widget_wavelengths.start) or (wavelength[-1] != widget_wavelengths.end) or ((wavelength[2]-wavelength[1]) != widget_wavelengths.step):
		for i,j in enumerate(range(widget_wavelengths.end - widget_wavelengths.start)):
			cds.wavelength.data['x'][
	'''

	''' input handlding '''
	if (widget_grating_types.active == 1):
		delta_lambda_default = 368
	else: # defaults to low resolution
		delta_lambda_default = 222

	delta_lambda = delta_lambda_default * widget_slit_width.value / 0.7

	if (0 in widget_plot_types.active) and (1 in widget_plot_types.active):
		channel = 'both'
	elif (0 in widget_plot_types.active):
		channel = 'blue'
	elif (1 in widget_plot_types.active):
		channel = 'red'
	else:
		channel = 'both'

	if (widget_moon_days.active == 0):
		sky_background = skyfile_00
	elif (widget_moon_days.active == 1):
		sky_background = skyfile_03
	elif (widget_moon_days.active == 2):
		sky_background = skyfile_05
	elif (widget_moon_days.active == 3):
		sky_background = skyfile_07
	elif (widget_moon_days.active == 4):
		sky_background = skyfile_10
	elif (widget_moon_days.active == 5):
		sky_background = skyfile_14
	else:
		raise ValueError('{} Invalid number of days since new moon ({})'.format(string_prefix,widget_moon_days.active))

	if (widget_telescope_sizes.active == 0):
		area = 368
	else:
		area = 222

	plot_step = wavelength[2] - wavelength[1]
	if ((wavelength[0] >= (3200-plot_step)) and (wavelength[-1] <= (10360+plot_step))): # picks up either unit
		plot_step = wavelength[2] - wavelength[1]
		#plot_step = delta_lambda / 3 default
	else:
		raise ValueError('{} Invalid wavelength range!'.format(string_prefix))

	# print('{} delta lambda: {} Angstrom, binned pixel scale {} Angstrom/px'.format(string_prefix,delta_lambda,plot_step)) # for debugging
	if (widget_object_types.active == 0):
		if (widget_types_types.active == 0):
			object_type = o5v
		elif (widget_types_types.active == 1):
			object_type = b0v
		elif (widget_types_types.active == 2):
			object_type = b57v
		elif (widget_types_types.active == 3):
			object_type = a0v
		elif (widget_types_types.active == 4):
			object_type = a5v
		elif (widget_types_types.active == 5):
			object_type = f0v
		elif (widget_types_types.active == 6):
			object_type = f5v
		elif (widget_types_types.active == 7):
			object_type = g0v
		elif (widget_types_types.active == 8):
			object_type = g5v
		elif (widget_types_types.active == 9):
			object_type = k0v
		elif (widget_types_types.active == 10):
			object_type = k5v
		elif (widget_types_types.active == 11):
			object_type = m0v
		elif (widget_types_types.active == 12):
			object_type = m5v
		else:
			pass

	elif (widget_object_types.active == 1):
		if (widget_galaxy_type.active == 0):
			object_type = sb1
		elif (widget_galaxy_type.active == 1):
			object_type = sb2
		elif (widget_galaxy_type.active == 2):
			object_type = sb3
		elif (widget_galaxy_type.active == 3):
			object_type = sb4
		elif (widget_galaxy_type.active == 4):
			object_type = sb5
		elif (widget_galaxy_type.active == 5):
			object_type = sb6
		elif (widget_galaxy_type.active == 6):
			object_type = s0
		elif (widget_galaxy_type.active == 7):
			object_type = sa
		elif (widget_galaxy_type.active == 8):
			object_type = sb
		elif (widget_galaxy_type.active == 9):
			object_type = sc
		elif (widget_galaxy_type.active == 10):
			object_type = bulge
		elif (widget_galaxy_type.active == 11):
			object_type = ellipticals
		elif (widget_galaxy_type.active == 12):
			object_type = lbg_all_flam
		else:
			pass
	else:
		object_type = a5v # default

	if (widget_filters.active == 0):
		selected_filter = ff_photon_ux
	elif (widget_filters.active == 1):
		selected_filter = ff_photon_b
	elif (widget_filters.active == 2):
		selected_filter = ff_photon_v
	elif (widget_filters.active == 3):
		selected_filter = ff_photon_r
	elif (widget_filters.active == 4):
		selected_filter = ff_photon_i
	elif (widget_filters.active == 5):
		selected_filter = ff_u
	elif (widget_filters.active == 6):
		selected_filter = ff_g
	elif (widget_filters.active == 7):
		selected_filter = ff_r
	elif (widget_filters.active == 8):
		selected_filter = ff_i
	elif (widget_filters.active == 9):
		selected_filter = ff_z
	else:
		selected_filter = ff_r # default

	if (widget_mag_type.active == 0): # vega
		mag_sys_opt = 'vega'
	else:
		mag_sys_opt = 'ab'


	filter_min = min(selected_filter.data['x'])
	filter_max = max(selected_filter.data['x'])

	lambda_min = filter_min
	lambda_max = filter_max

	plot_step = wavelength[2] - wavelength[1]
	lambda_A = range(lambda_min,lambda_max,plot_step)

	object_x = object_type[0] * (1 + widget_redshift.value)
	object_y = object_type[1]

	cds_spectres.data['n'] = lambda_A
	cds_spectres.data['o'] = object_x
	cds_spectres.data['s'] = object_y
	cds_spectres.change.emit()
	flux_A = cds_spectres.data['f']
	#flux_A = spectres(lambda_A,object_x,object_y)
	print(flux_A[5])
	lambda_A[0] = lambda_A[0] + plot_step
	lambda_A[-1] = lambda_A[-1] - plot_step

	ftrans = interpolate.interp1d(selected_filter[0],selected_filter[1], kind='cubic')
	trans = ftrans(lambda_A) #spectres(lambda_A,selected_filter[0],selected_filter[1])

	atmo_ext_x = atmo_ext_file.data['x']
	atmo_ext_y = atmo_ext_file.data['y']

	extinction = spectres(lambda_A,atmo_ext_x,atmo_ext_y)

	flux = flux_A * 1e10
	_lambda = lambda_A / 1e10

	num_zeros = 0
	for lux in flux:
		if (lux is None) or (lux is 0):
			num_zeros += 1
			lux = 0

	if (num_zeros >= (flux.shape[0]/5)):
		if (num_zeros == flux.shape[0]):
			print('No flux in this bandpass!')
			output_flux = [0 for i in range(wavelength.shape[0])] # np.zeros
		else:
			percent_zeros = (num_zeros / flux.shape[0]) * 100
			print('{}% of this bandpass has zero flux'.format(percent_zeros))

	if (mag_sys_opt == 'vega'):
		flux_vega = spectres(wavelength,vega_file.data['x'],vega_file.data['x']) * 1e10
		print(vega[0])
		mag_model = -2.5 * [math.log10(i) for i in (math.fsum(flux * extinction * _lambda * trans) / math.fsum(flux_vega * trans * _lambda * extinction))] + 0.03
	elif (mag_sys_opt == 'ab'):
		#mag_model = -48.60 - 2.5 * np.log10(np.divide(math.fsum(np.multiply(flux,trans).multiply(extinction).multiply(_lambda)),math.fsum(trans.multiply(_lambda) * extinction[1]).multiply(const.c.value/(np.square(_lambda)))))
		mag_model = -48.6 - 2.5 * [math.log10(i) for i in (math.fsum(flux * trans * extinction *_lambda) / math.fsum(trans * _lambda * extinction * (astroconst.c.value/[i**2 for i in (_lambda)])))]
	else:
		raise ValueError('Invalid magnitude system option!')

	del_mag = mag - mag_model
	star_x = object_x
	print(type(del_mag))
	star_y = np.multiply(object_y,10 ** -del_mag/2.5)

	old_res = star_x[2] - star_x[1]
	if (old_res < plot_step):
		flux_y = spectres(wavelength,star_x,(star_y*1e-03)) # ergs s-1 cm-2 A-1 to J s-1 m-2 A-1
	else:
		flux_y = spectres(wavelength,star_x,(star_y*1e-03))

	flux = flux_y
	power = flux * area * widget_exposure_time.value * plot_step
	counts = np.divide(np.divide(power,np.divide((astroconst.h.value * astroconst.c.value),wavelength)),1e10)

	''' subtract light lost to components '''
	# seeing
	_sigma = widget_seeing.value / gaussian_sigma_to_fwhm
	funx = lambda x: (1/(_sigma * math.sqrt(2 * math.pi)))*math.exp((-(x ** 2)/(2 * (_sigma ** 2))))
	percent_u,percent_err_u = integrate.quad(funx,(-widget_slit_width.value/2),(widget_slit_width.value/2))
	percent_l,percent_err_l = integrate.quad(funx,(-widget_seeing.value/2),(widget_seeing.value/2))
	percent = percent_u * percent_l # can use error if you add it later...
	extension = widget_seeing.value * widget_slit_width.value

	# sky background
	sky_x = sky_background.data['x'] * 1e4
	sky_y = sky_background.data['y'] / 1e4
	old_res = sky_x[2] - sky_x[1]
	_sigma = delta_lambda / gaussian_sigma_to_fwhm
	_x = range((-5*_sigma),(5*_sigma),old_res)
	degrade = funx(_x)/np.trapz(funx(_x))
	sky_y = convolve_fft(sky_y,degrade)
	sky_flux = spectres(wavelength,sky_x,sky_y)
	counts_noise = np.multiply(np.multiply(sky_flux,extension),(area * widget_exposure_time.value * plot_step))

	# dichroic
	if (channel is 'blue') or (channel is 'both'):
		fblue_dichro = interpolate.interp1d(dichroic_file.data['x'],dichroic_file.data['y1'], kind='cubic')
		blue_dichro = fblue_dichro(wavelength)
	if (channel is 'red') or (channel is 'both'):
		fred_dichro = interpolate.interp1d(dichroic_file.data['x'],dichroic_file.data['y2'], kind='cubic')
		red_dichro = fred_dichro(wavelength)

	# grating
	if (channel is 'blue') or (channel is 'both'):
		blue_grating = spectres(wavelength,(grating_blue.data['x']*10),grating_blue.data['y'])
	if (channel is 'red') or (channel is 'both'):
		red_grating = spectres(wavelength,(grating_red.data['x']*10),grating_red.data['y'])

	# ccd
	if (channel is 'blue') or (channel is 'both'):
		fblue_ccd = interpolate.interp1d((ccd_blue.data['x']*10),ccd_blue.data['y'], kind='cubic')
		blue_ccd = fblue_ccd(wavelength)
	if (channel is 'red') or (channel is 'both'):
		fred_ccd = interpolate.interp1d((ccd_red.data['x']*10),ccd_red.data['y'], kind='cubic')
		red_ccd = fred_ccd(wavelength)

	# mirror bouncing, apply twice
	fmirror = interpolate.interp1d(mirror_file.data['x'],mirror_file.data['y'], kind='cubic')
	mirror = fmirror(wavelength)

	# read noise
	spectral_resolution = math.ceil((widget_slit_width.value/(0.7/12))/2)*2 #px (ceil()/2)*2 to round up to next even integer
	spatial_resolution = math.ceil((widget_seeing.value/(0.7/12))/2)*2 #px (ceil()/2)*2 to round up to next even integer 
	extent = widget_seeing.value * widget_slit_width.value
	npix = extent/(0.7/12)**2
	
	bin_size = widget_binned_pixel_scale.active + 1 # default 2x2 binning

	rn = 2 # default e-
	if (bin_size > 0) and (bin_size < 5):
		print('[ info ] : Pixel binning: ({}x{})'.format(bin_size,bin_size))
		readnoise = math.ceil(rn * spectral_resolution * spatial_resolution / (bin_size**2))
		print('[ info ] : Extent: {} arcsec^2\n[ info ] : num pixels/resel: {} px\n[ info ] : spectral resolution: {} px\n[ info ] : spatial resolution: {} px'.format(extent,int(math.ceil(npix)),spectral_resolution,spatial_resolution))
	else:
		raise ValueError('{} Invalid pixel binning option ({})'.format(string_prefix,bin_size))

	extinction = spectres(wavelength,atmo_ext_file.data['x'],atmo_ext_file.data['y']) # since not passed from function, just use global in real version

	''' calculations '''

	# signal
	if (channel == 'blue') or (channel == 'both'):
		blue_total_eff = np.multiply(np.multiply(dichroic_file.data['y1'],grating_blue.data['y']),np.multiply((ccd_blue.data['y'] * (coating_eff_blue * extinction)),np.square(mirror_file.data['y'])))
		blue_signal = np.multiply((counts * percent), blue_total_eff)
	if (channel == 'red') or (channel == 'both'):
		red_total_eff = np.multiply(np.multiply(dichroic_file.data['y2'],grating_red.data['y']),np.multiply((ccd_red.data['y'] * (coating_eff_red * extinction)),np.square(mirror_file.data['y'])))
		red_signal = np.multiply((counts * percent), red_total_eff)

	# noise
	if (channel == 'blue') or (channel == 'both'):
		blue_total_eff_noise = np.multiply(np.multiply(dichroic_file.data['y1'],grating_blue.data['y']),(ccd_blue.data['y'] * np.square(mirror.file['y']) * coating_eff_blue))
		blue_noise = np.multiply(counts_noise,blue_total_eff_noise)
	if (channel == 'red') or (channel == 'both'):
		red_total_eff_noise = np.multiply(np.multiply(dichroic_file.data['y1'],grating_red.data['y']),(ccd_red.data['y'] * np.square(mirror.file['y']) * coating_eff_red))
		red_noise = np.multiply(counts_noise,red_total_eff_noise)

	if (channel == 'blue') or (channel == 'both'):
		snr_blue = np.divide(blue_signal,np.sqrt(blue_signal + blue_noise + np.square(readnoise)))
	if (channel == 'red') or (channel == 'both'):
		snr_red = np.divide(red_signal,np.sqrt(red_signal + red_noise + np.square(readnoise)))

	if (channel == 'blue') or (channel == 'both'):
		sigma_blue = np.sqrt(blue_signal + blue_noise + np.square(readnoise))
	if (channel == 'red') or (channel == 'both'):
		sigma_red = np.sqrt(red_signal + red_noise + np.square(readnoise))

	if (channel == 'blue') or (channel == 'both'):
		error_blue = np.random.normal(loc=0, scale=sigma_blue,size=len(snr_blue))
	if (channel == 'red') or (channel == 'both'):
		error_red = np.random.normal(loc=0, scale=sigma_red,size=len(snr_red))


	''' pre-plotting '''
	plot_y_blue,plot_y_red = [],[]
	if (tabs.active == 0):
		title = 'Signal-to-Noise Ratio'
		labels = ['Angstrom','Signal-to-Noise Ratio']
		if (channel == 'blue') or (channel == 'both'):
			plot_y_blue = snr_blue
		if (channel == 'red') or (channel == 'both'):
			plot_y_red = snr_red
	elif (tabs.active == 1):
		labels = ['Angstrom','Counts (#) per pixel']
		if widget_withnoise.active == 1:
			title = 'Observed Spectrum (with noise)'
			if (channel == 'blue') or (channel == 'both'):
				plot_y_blue = np.add(blue_signal,error_blue)
			if (channel == 'red') or (channel == 'both'):
				plot_y_red = np.add(red_signal,error_red)
		else:
			title = 'Observed Spectrum (without noise)'
			if (channel == 'blue') or (channel == 'both'):
				plot_y_blue = blue_signal
			if (channel == 'red') or (channel == 'both'):
				plot_y_red = red_signal
	elif (tabs.active == 2):
		title = 'Observed Sky Background'
		labels = ['Angstrom','Counts (#) per pixel']
		if (channel == 'blue') or (channel == 'both'):
			plot_y_blue = blue_noise
		if (channel == 'red') or (channel == 'both'):
			plot_y_red = red_noise
	elif (tabs.active == 3):
		title = 'Dichroic Throughput'
		labels = ['Angstrom','Throughput']
		if (channel == 'blue') or (channel == 'both'):
			plot_y_blue = blue_dichro
		if (channel == 'red') or (channel == 'both'):
			plot_y_red = red_dichro
	elif (tabs.active == 4):
		title = 'Grating Throughput'
		labels = ['Angstrom','Throughput']
		if (channel == 'blue') or (channel == 'both'):
			plot_y_blue = blue_grating
		if (channel == 'red') or (channel == 'both'):
			plot_y_red = red_grating
	elif (tabs.active == 5):
		title = 'CCD Quantum Efficiency'
		labels = ['Angstrom','QE']
		if (channel == 'blue') or (channel == 'both'):
			plot_y_blue = blue_ccd
		if (channel == 'red') or (channel == 'both'):
			plot_y_red = red_ccd
	elif (tabs.active == 6):
		title = 'Atmospheric Extinction'
		labels = ['Angstrom','Throughput']
		plot_y_blue = spectres(wavelength,atmo_ext_x,atmo_ext_y)
		plot_y_red = plot_y_blue


		cds_blue.data['xr'] = wavelength
		cds_red.data['xr'] = wavelength

		cds_blue.data['yb'] = plot_y_blue
		cds_red.data['yr'] = plot_y_red
		cds_wavelength.data['x'] = wavelength

		cds_blue.change.emit()
		cds_red.change.emit()


# glyphs
gly_red = glyphs.Line(x="xr",y="yr",line_color='red')
gly_blue = glyphs.Line(x="xb",y="yb",line_color='blue')
p0.add_glyph(cds_blue,gly_blue)
p0.add_glyph(cds_red,gly_red)
p1.add_glyph(cds_blue,gly_blue)
p1.add_glyph(cds_red,gly_red)
p2.add_glyph(cds_blue,gly_blue)
p2.add_glyph(cds_red,gly_red)
p3.add_glyph(cds_blue,gly_blue)
p3.add_glyph(cds_red,gly_red)
p4.add_glyph(cds_blue,gly_blue)
p4.add_glyph(cds_red,gly_red)
p5.add_glyph(cds_blue,gly_blue)
p5.add_glyph(cds_red,gly_red)
p6.add_glyph(cds_blue,gly_blue)
p6.add_glyph(cds_red,gly_red)


def plot_types_callback(gly_blue=gly_blue,gly_red=gly_red):
	if (0 in cb_obj.active):
		gly_blue.line_alpha = 0.5
	elif (0 not in cb_obj.active):
		gly_blue.line_alpha = 0.0
	if (1 in cb_obj.active):
		gly_red.line_alpha = 0.5
	elif (1 not in cb_obj.active):
		gly_red.line_alpha = 0.0


spectres_callback = CustomJS.from_py_func(spectres)
#cds_spectres.js_on_change('data',spectres_callback)
cds_spectres.callback = spectres_callback

gen_callback = CustomJS.from_py_func(general_callback)
widget_object_types.callback = gen_callback
widget_types_types.callback = gen_callback
widget_galaxy_type.callback = gen_callback
widget_grating_types.callback = gen_callback
widget_binned_pixel_scale.callback = gen_callback
widget_wavelengths.callback = gen_callback
widget_plot_types.callback = CustomJS.from_py_func(plot_types_callback)


# final panel building
widget_group_one = widgetbox(children=[widget_telescope_sizes,widget_object_types,widget_types_types,widget_galaxy_type])
widget_group_two = layout([[widget_mag_input],[widget_filters,widget_mag_type]])
widget_group_three = widgetbox(children=[widget_grating_types,widget_redshift,widget_exposure_time,widget_seeing,widget_slit_width,widget_moon_days_header,widget_moon_days,widget_wavelengths,widget_binned_pixel_scale,widget_plot_types]) # removed widget_plot_step and widget_binned_pixel_scale_mode
widgets = column(children=[widget_group_one,widget_group_two,widget_group_three],width=dfs.default_toolbar_width)
inputs = row(children=[widgets,tabs],sizing_mode='scale_height')

l = layout([[widget_header],[inputs]])
show(l)

print('Build completed in {} seconds.'.format(time.time()-time_start))