#!/usr/bin/python3
'''
	GMACS ETC: UIC
		(user-interface creator)
			builds the etc webpage

'''
import os
import time
import json
import math
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy import units as u
from util import config as cfg
from util import defaults as dfs
from bokeh.plotting import figure
from bokeh.embed import autoload_static
from bokeh.io import output_file,show,save
from bokeh.models.callbacks import CustomJS
from bokeh.layouts import widgetbox, column, row, layout
from bokeh.models import ColumnDataSource, glyphs, ranges, HoverTool, CustomJSHover
from bokeh.models.widgets import Dropdown, RadioButtonGroup, CheckboxButtonGroup, Slider, TextInput, RangeSlider, Paragraph, Panel, Tabs, Div

output_file(dfs.string_pagename, title=dfs.site_title, mode='cdn')

'''
	GMACS ETC: ODC
		omni-dimensional creator
			assuage my confusion about creation
				and the universe

'''
print('Building...')
time_start = time.time()

# explicit... not very elegant, but pandas don't like loops
galaxy_sb1 = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[0]),sep='\s+',skiprows=1))
galaxy_sb2 = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[1]),sep='\s+',skiprows=1))
galaxy_sb3 = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[2]),sep='\s+',skiprows=1))
galaxy_sb4 = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[3]),sep='\s+',skiprows=1))
galaxy_sb5 = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[4]),sep='\s+',skiprows=1))
galaxy_sb6 = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[5]),sep='\s+',skiprows=1))
galaxy_s0 = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[6]),sep='\s+',skiprows=1))
galaxy_sa = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[7]),sep='\s+',skiprows=1))
galaxy_sb = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[8]),sep='\s+',skiprows=1))
galaxy_sc = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[9]),sep='\s+',skiprows=1))
galaxy_bulge = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[10]),sep='\s+',skiprows=1))
galaxy_ellipticals = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[11]),sep='\s+',skiprows=1))
galaxy_lbg_all_flam = ColumnDataSource(pd.read_csv(os.path.join(cfg.galaxy_types_path,cfg.galaxy_types_files[12]),sep='\s+',skiprows=1))

panda_05v = pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[0]),sep='\s+',skiprows=1)
star_o5v = ColumnDataSource(panda_05v) # trust
star_b0v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[1]),sep='\s+',skiprows=1))
star_b57v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[2]),sep='\s+',skiprows=1))
star_a0v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[3]),sep='\s+',skiprows=1))
star_a5v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[4]),sep='\s+',skiprows=1))
star_f0v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[5]),sep='\s+',skiprows=1))
star_g0v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[6]),sep='\s+',skiprows=1))
star_g5v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[7]),sep='\s+',skiprows=1))
star_k0v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[8]),sep='\s+',skiprows=1))
star_k5v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[9]),sep='\s+',skiprows=1))
star_m0v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[10]),sep='\s+',skiprows=1))
star_m5v = ColumnDataSource(pd.read_csv(os.path.join(cfg.star_types_path,cfg.star_types_files[11]),sep='\s+',skiprows=1))

filter_photonux = ColumnDataSource(pd.read_csv(os.path.join(cfg.filter_path,cfg.filter_files[0]),sep='\s+',skiprows=1))
filter_photonb = ColumnDataSource(pd.read_csv(os.path.join(cfg.filter_path,cfg.filter_files[1]),sep='\s+',skiprows=1))
filter_photonv = ColumnDataSource(pd.read_csv(os.path.join(cfg.filter_path,cfg.filter_files[2]),sep='\s+',skiprows=1))
filter_photonr = ColumnDataSource(pd.read_csv(os.path.join(cfg.filter_path,cfg.filter_files[3]),sep='\s+',skiprows=1))
filter_photoni = ColumnDataSource(pd.read_csv(os.path.join(cfg.filter_path,cfg.filter_files[4]),sep='\s+',skiprows=1))
filter_u = ColumnDataSource(pd.read_csv(os.path.join(cfg.filter_path,cfg.filter_files[5]),sep='\s+',skiprows=1))
filter_g = ColumnDataSource(pd.read_csv(os.path.join(cfg.filter_path,cfg.filter_files[6]),sep='\s+',skiprows=1))
filter_r = ColumnDataSource(pd.read_csv(os.path.join(cfg.filter_path,cfg.filter_files[7]),sep='\s+',skiprows=1))
filter_i = ColumnDataSource(pd.read_csv(os.path.join(cfg.filter_path,cfg.filter_files[8]),sep='\s+',skiprows=1))
filter_z = ColumnDataSource(pd.read_csv(os.path.join(cfg.filter_path,cfg.filter_files[9]),sep='\s+',skiprows=1))

skyfile_00d = ColumnDataSource(pd.read_csv(os.path.join(cfg.skyfiles_path,cfg.skyfiles_files[0]),sep='\s+',skiprows=1))
skyfile_03d = ColumnDataSource(pd.read_csv(os.path.join(cfg.skyfiles_path,cfg.skyfiles_files[1]),sep='\s+',skiprows=1))
skyfile_07d = ColumnDataSource(pd.read_csv(os.path.join(cfg.skyfiles_path,cfg.skyfiles_files[2]),sep='\s+',skiprows=1))
skyfile_10d = ColumnDataSource(pd.read_csv(os.path.join(cfg.skyfiles_path,cfg.skyfiles_files[3]),sep='\s+',skiprows=1))
panda_14d = pd.read_csv(os.path.join(cfg.skyfiles_path,cfg.skyfiles_files[4]),sep='\s+',skiprows=1)
skyfile_14d = ColumnDataSource(panda_14d)

efficiency_grating_red = pd.read_csv(os.path.join(cfg.efficiency_path,cfg.efficiency_grating_files[0]),sep=',',skiprows=1)
efficiency_grating_blue = pd.read_csv(os.path.join(cfg.efficiency_path,cfg.efficiency_grating_files[1]),sep=',',skiprows=1)
efficiency_ccd_red = pd.read_csv(os.path.join(cfg.efficiency_path,cfg.efficiency_ccd_files[0]),sep='\s+',skiprows=1)
efficiency_ccd_blue = pd.read_csv(os.path.join(cfg.efficiency_path,cfg.efficiency_ccd_files[0]),sep='\s+',skiprows=1)

dichroic = pd.read_csv(os.path.join(cfg.dichroic_path,cfg.dichroic_files[0]),sep='\s+',skiprows=1)
atmo_ext = pd.read_csv(os.path.join(cfg.atmo_ext_path,cfg.atmo_ext_files[0]),sep=',',skiprows=1)

wavelength = np.arange(start=dfs.default_start_wavelength,stop=dfs.default_stop_wavelength,step=dfs.default_plot_step_step)
cds_wavelength = ColumnDataSource(dict(wavelength=wavelength))

# efficiency dots
#print(efficiency_grating_red[str(efficiency_grating_red.columns[0])])
#TODO arg
#grating_red_1 = np.dot(efficiency_grating_red[str(efficiency_grating_red.columns[0])],10)
grating_red_1 = efficiency_grating_red[str(efficiency_grating_red.columns[0])]
grating_red_2 = efficiency_grating_red[str(efficiency_grating_red.columns[1])]
#grating_blue_1 = np.dot(efficiency_grating_blue[str(efficiency_grating_blue.columns[0])],10)
grating_blue_1 = efficiency_grating_blue[str(efficiency_grating_blue.columns[0])]
grating_blue_2 = efficiency_grating_blue[str(efficiency_grating_blue.columns[1])]

#ccd_efficiency_red_1 = np.dot(efficiency_ccd_red[0],10)
ccd_efficiency_red_1 = efficiency_ccd_red[str(efficiency_ccd_red.columns[0])]
ccd_efficiency_red_2 = efficiency_ccd_red[str(efficiency_ccd_red.columns[1])]
#ccd_efficiency_blue_1 = np.dot(efficiency_ccd_blue[0],10)
ccd_efficiency_blue_1 = efficiency_ccd_blue[str(efficiency_ccd_blue.columns[0])]
ccd_efficiency_blue_2 = efficiency_ccd_blue[str(efficiency_ccd_blue.columns[1])]

#dichro_x = np.dot(dichroic[0],10) # wavelength in Angstroms
dichro_x = dichroic[str(dichroic.columns[0])]
dichro_y1 = dichroic[str(dichroic.columns[1])] # reflectivity, blue channel
dichro_y2 = dichroic[str(dichroic.columns[2])]

atmo_ext_x = atmo_ext[str(atmo_ext.columns[0])]
atmo_ext_y = atmo_ext[str(atmo_ext.columns[1])]


'''
# divide (n,x,y) into {(n,x),(n,y)}
coalesced_object_x = {}
coalesced_object_y = {}
for i,coalesced_galaxy_type in tqdm(enumerate(coalesced_galaxy_types),desc='Shaping data',ncols=0):
	coalesced_object_x[i] = coalesced_galaxy_types[i].a
	coalesced_object_y[i] = coalesced_galaxy_types[i].b
'''

#pandas_dataframe = ColumnDataSource().add(coalesced_galaxy_types,coalesced_star_types)
#pandas_dataframe = ColumnDataSource()
#pandas_dataframe = pandas_dataframe.add(data=[frame for frame in frames])
#pandas_dataframe = ColumnDataSource(pandas_dataframe)


'''
	ODC User Interface
		manages all interactivity

'''
# radio button groups
widget_telescope_sizes = RadioButtonGroup(labels=dfs.string_telescope_sizes, active=0, name=dfs.string_widget_names[0])
widget_object_types = RadioButtonGroup(labels=dfs.string_object_types, active=dfs.default_object_types, name=dfs.string_widget_names[1])
widget_galaxy_type = Dropdown(label=dfs.string_object_types_types[1],menu=dfs.string_galaxy_types, name=dfs.string_widget_names[2])
widget_galaxy_type.disabled = True # change if galaxies default
widget_mag_type = RadioButtonGroup(labels=dfs.string_magnitude_three,active=0, name=dfs.string_widget_names[3])
widget_grating_types = RadioButtonGroup(labels=dfs.string_grating_types, active=0, name=dfs.string_widget_names[4])
widget_moon_days = RadioButtonGroup(labels=dfs.string_moon_days, active=0, name=dfs.string_widget_names[5]) # TODO: need title for this, would be string_title[7]
widget_binned_pixel_scale_mode = RadioButtonGroup(labels=dfs.string_binned_pixel_scale_modes, active=0)
widget_binned_pixel_scale = RadioButtonGroup(name=dfs.string_binned_pixel_scale_header, labels=dfs.string_binned_pixel_scale_labels,active=0)
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


'''
plot it
- figures are the... yeah if you're reading this, then you know what that means
- tabN 'Panels' are panel spaces that are flipped between
- the 'tabs' variable builds the tabs bar and points them toward their respective panels

'''
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


'''
	Python functions for conversion to JavaScript
		Remember that none of these functions are run on the page;
			These functions are only here to be converted to JavaScript,
				Which are then linked to the widgets as JavaScript callbacks.

				omggg print(star_o5v.data[list(star_o5v.data)[0]])
'''
"""
SpectRes: A fast spectral new_spec_wavs function.
Copyright (C) 2017  A. C. Carnall
GNU General Public License
A couple custom mods for this application
"""
def spectres(new_spec_wavs=cds_wavelength.data['wavelength'], old_spec_wavs=cds_snr_red.data['xr'], spec_fluxes=cds_snr_red.data['yr'], spec_errs=None):
    spec_lhs = np.zeros(np.asarray(old_spec_wavs).shape[0])
    spec_widths = np.zeros(np.asarray(old_spec_wavs).shape[0])
    spec_lhs = np.zeros(np.asarray(old_spec_wavs).shape[0])
    spec_lhs[0] = old_spec_wavs[0] - (old_spec_wavs[1] - old_spec_wavs[0])/2
    spec_widths[-1] = (old_spec_wavs[-1] - old_spec_wavs[-2])
    spec_lhs[1:] = (np.asarray(old_spec_wavs)[1:] + np.asarray(old_spec_wavs)[:-1])/2
    spec_widths[:-1] = spec_lhs[1:] - spec_lhs[:-1]
    filter_lhs = np.zeros(new_spec_wavs.shape[0]+1)
    filter_widths = np.zeros(new_spec_wavs.shape[0])
    filter_lhs[0] = new_spec_wavs[0] - (new_spec_wavs[1] - new_spec_wavs[0])/2
    filter_widths[-1] = (new_spec_wavs[-1] - new_spec_wavs[-2])
    filter_lhs[-1] = new_spec_wavs[-1] + (new_spec_wavs[-1] - new_spec_wavs[-2])/2
    filter_lhs[1:-1] = (new_spec_wavs[1:] + new_spec_wavs[:-1])/2
    filter_widths[:-1] = filter_lhs[1:-1] - filter_lhs[:-2]
    if filter_lhs[0] < spec_lhs[0] or filter_lhs[-1] > spec_lhs[-1]:
        raise ValueError("spectres: The new wavelengths specified must fall within the range of the old wavelength values.")
 #b =  array.slice().map( function(row){ return row.slice(); });
    resampled_fluxes = np.zeros(spec_fluxes[0].shape + new_spec_wavs.shape)
    if spec_errs is not None:
        if spec_errs.shape != spec_fluxes.shape:
            raise ValueError("If specified, spec_errs must be the same shape as spec_fluxes.")
        else:
            resampled_fluxes_errs = np.copy(resampled_fluxes)
    start = 0
    stop = 0
    for j in range(new_spec_wavs.shape[0]):
        while spec_lhs[start+1] <= filter_lhs[j]:
            start += 1
        while spec_lhs[stop+1] < filter_lhs[j+1]:
            stop += 1
        if stop == start:
            resampled_fluxes[j] = spec_fluxes[start]
            if spec_errs is not None:
                resampled_fluxes_errs[j] = spec_errs[start]
        else:
            start_factor = (spec_lhs[start+1] - filter_lhs[j])/(spec_lhs[start+1] - spec_lhs[start])
            end_factor = (filter_lhs[j+1] - spec_lhs[stop])/(spec_lhs[stop+1] - spec_lhs[stop])
            spec_widths[start] *= start_factor
            spec_widths[stop] *= end_factor
            resampled_fluxes[j] = np.sum(spec_widths[start:stop+1]*spec_fluxes[start:stop+1], axis=-1)/np.sum(spec_widths[start:stop+1])
            if spec_errs is not None:
                resampled_fluxes_errs[j] = np.sqrt(np.sum((spec_widths[start:stop+1]*spec_errs[start:stop+1])**2, axis=-1))/np.sum(spec_widths[start:stop+1])
            spec_widths[start] /= start_factor
            spec_widths[stop] /= end_factor
    if spec_errs is not None:
        return resampled_fluxes, resampled_fluxes_errs
    else: 
        return resampled_fluxes

spect_res = CustomJS.from_py_func(spectres)
def fun_callback(galaxy_sb1=galaxy_sb1,galaxy_sb2=galaxy_sb2,galaxy_sb3=galaxy_sb3,galaxy_sb4=galaxy_sb4,galaxy_sb5=galaxy_sb5,galaxy_sb6=galaxy_sb6,galaxy_s0=galaxy_s0,galaxy_sa=galaxy_sa,
				galaxy_sb=galaxy_sb,galaxy_sc=galaxy_sc,galaxy_bulge=galaxy_bulge,galaxy_ellipticals=galaxy_ellipticals,galaxy_lbg_all_flam=galaxy_lbg_all_flam,
				star_o5v=star_o5v,star_b0v=star_b0v,star_b57v=star_b57v,star_a0v=star_a0v,star_a5v=star_a5v,star_f0v=star_f0v,
				star_g0v=star_g0v,star_g5v=star_g5v,star_k0v=star_k0v,star_k5v=star_k5v,star_m0v=star_m0v,star_m5v=star_m5v,
				filter_photonux=filter_photonux,filter_photonb=filter_photonb,filter_photonv=filter_photonv,filter_photonr=filter_photonr,filter_photoni=filter_photoni,
				filter_u=filter_u,filter_g=filter_g,filter_r=filter_r,filter_i=filter_i,filter_z=filter_z,
				skyfile_00d=skyfile_00d,skyfile_03d=skyfile_03d,skyfile_07d=skyfile_07d,skyfile_10d=skyfile_10d,skyfile_14d=skyfile_14d,
				cds_snr_red=cds_snr_red,cds_snr_blue=cds_snr_blue,cds_os_noise_red=cds_os_noise_red,cds_os_noise_blue=cds_os_noise_blue,
				cds_os_nonoise_red=cds_os_noise_red,cds_os_nonoise_blue=cds_os_nonoise_blue,cds_wavelength=cds_wavelength,
				spect_res=spect_res,
				tabs=tabs,widget_object_types=widget_object_types,widget_galaxy_type=widget_galaxy_type,widget_types_types=widget_types_types,widget_mag_type=widget_mag_type,
				widget_wavelengths=widget_wavelengths):
	print('[ETC] Call from ' + cb_obj.name) # for debugging
	if (tabs.active == 0): # signal-to-noise

		if (cb_obj.name == widget_object_types.name):
			# handling stellar classification changes
			if (widget_object_types.active == 0):
				widget_types_types.disabled = False
				widget_galaxy_type.disabled = True
				widget_mag_type.disabled = False
				print('[ETC] Object type: Stars')
				# if something selected previously...
				if (widget_types_types.value != None):
					print('[ETC] Stellar classification: {}'.format(widget_types_types.value))

			# handling galactic classification changes
			else:
				widget_types_types.disabled = True
				widget_galaxy_type.disabled = False
				widget_mag_type.disabled = True
				print('[ETC] Object type: Galaxies')
				# if something selected previously...
				if (widget_galaxy_type.value != None):
					print('[ETC] Galactic classification: {}'.format(widget_galaxy_type.value))

		# a freshly selected galatic classification
		elif (cb_obj.name == widget_galaxy_type.name):
			widget_galaxy_type.label = widget_galaxy_type.menu[widget_galaxy_type.value][0] # trust
			# explicit TODO not explicit
			if (widget_galaxy_type.value[1] == 0):
				_gal_file = galaxy_sb1
			elif (widget_galaxy_type.value[1] == 1):
				_gal_file = galaxy_sb2
			elif (widget_galaxy_type.value[1] == 2):
				_gal_file = galaxy_sb3
			elif (widget_galaxy_type.value[1] == 3):
				_gal_file = galaxy_sb4
			elif (widget_galaxy_type.value[1] == 4):
				_gal_file = galaxy_sb5
			elif (widget_galaxy_type.value[1] == 5):
				_gal_file = galaxy_sb6
			elif (widget_galaxy_type.value[1] == 6):
				_gal_file = galaxy_s0
			elif (widget_galaxy_type.value[1] == 7):
				_gal_file = galaxy_sa
			elif (widget_galaxy_type.value[1] == 8):
				_gal_file = galaxy_sb
			elif (widget_galaxy_type.value[1] == 9):
				_gal_file = galaxy_sc
			elif (widget_galaxy_type.value[1] == 10):
				_gal_file = galaxy_bulge
			elif (widget_galaxy_type.value[1] == 11):
				_gal_file = galaxy_ellipticals
			elif (widget_galaxy_type.value[1] == 11):
				_gal_file = galaxy_lbg_all_flam

				#TODO: replace 0.1, use variable global_plot_step or something
			_wavelength = [widget_wavelengths.value[0]]+[(widget_wavelengths.value[0]+(0.1*i)) for i in range(int((widget_wavelengths.value[1]-widget_wavelengths.value[0])/0.1))]+[widget_wavelengths.value[1]]
			
			(cds_snr_red.data)['xr'],(cds_snr_red.data)['yr'] = _gal_file.data[list(_gal_file.data)[0]],_gal_file.data[list(_gal_file.data)[1]]
			
			#cds_snr_red

			print('[ETC] Galaxy type: {}'.format(widget_galaxy_type.label))

		# a freshly selected stellar classification
		elif (cb_obj.name == widget_types_types.name):
			widget_types_types.label = widget_types_types.value
			print('[ETC] Star type: {}'.format(widget_types_types.label))


	elif (tabs.active == 1): # spec w/o noise
		pass
	elif (tabs.active == 2): # spec w/ noise
		pass
	elif (tabs.active == 3): # sky background
		pass
	else:
		print('[ETC] Switched to predefined plot')

def plot_types_callback(gly_snr_red=gly_snr_red,gly_os_noise_red=gly_os_noise_red,gly_os_nonoise_red=gly_os_nonoise_red,
					gly_sky_red=gly_sky_red,gly_dichroic_red=gly_dichroic_red,gly_grating_red=gly_grating_red,gly_ccd_red=gly_ccd_red,
					gly_snr_blue=gly_snr_blue,gly_os_noise_blue=gly_os_noise_blue,gly_os_nonoise_blue=gly_os_nonoise_blue,
					gly_sky_blue=gly_sky_blue,gly_dichroic_blue=gly_dichroic_blue,gly_grating_blue=gly_grating_blue,gly_ccd_blue=gly_ccd_blue,tabs=tabs):
	red_glyphs = [gly_snr_red,gly_os_noise_red,gly_os_nonoise_red,gly_sky_red,gly_dichroic_red,gly_grating_red,gly_ccd_red]
	blue_glyphs = [gly_snr_blue,gly_os_noise_blue,gly_os_nonoise_blue,gly_sky_blue,gly_dichroic_blue,gly_grating_blue,gly_ccd_blue]
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

def read_noise(widget_seeing=widget_seeing,widget_slit_width=widget_slit_width,cds_noise=cds_noise,widget_binned_pixel_scale=widget_binned_pixel_scale,widget_grating_types=widget_grating_types):
	slit_size = widget_slit_width.value
	seeing = widget_seeing.value
	bps = widget_binned_pixel_scale.active + 1

	if (widget_grating_types.active == 1):
		delta_lambda_default = 1.4 # high res
	else:
		delta_lambda_default = 3.73 # low res, default
	if (bps == 0):
		bp_const = 12
	elif (bps == 2):
		bp_const = 4
	elif (bps == 3):
		bp_const = 3
	else:
		bp_const = 6
	rn = delta_lambda_default / bp_const
	spec_resl = int(slit_size/(0.7/12))
	spat_resl = int(seeing/(0.7/12))
	extent = seeing*slit_size
	npix = extent/(0.7/12)**2
	print('extent: ', extent, 'arcsec^2\t\t\t', 'num pixels: ',npix, 'px')
	print('spectral resolution: ', spec_resl, 'px\t\t', 'spatial resolution: ', spat_resl, 'px')
	readnoise = lambda binsize: int(rn * spec_resl * spat_resl / (binsize * binsize)) # e-
	_rn = readnoise(bps)
	print('bin {}x: {} e-'.format(bps,_rn))

# linkages
coalesced_callback = CustomJS.from_py_func(fun_callback) # only convert to JS once!
tabs.callback = coalesced_callback
widget_object_types.callback = coalesced_callback
widget_types_types.callback = coalesced_callback
widget_galaxy_type.callback = coalesced_callback
widget_plot_types.callback = CustomJS.from_py_func(plot_types_callback) # this can go straight in (unlike coalesced) since only one idget calls it; it only gets instanced once
eh = CustomJS.from_py_func(read_noise)
widget_grating_types.callback = eh
widget_binned_pixel_scale.callback = eh

# final panel building
widget_group_one = widgetbox(children=[widget_telescope_sizes,widget_object_types,widget_types_types,widget_galaxy_type])
widget_group_two = layout([[widget_mag_input],[widget_filters,widget_mag_type]])
widget_group_three = widgetbox(children=[widget_grating_types,widget_redshift,widget_exposure_time,widget_seeing,widget_slit_width,widget_moon_days_header,widget_moon_days,widget_wavelengths,widget_binned_pixel_scale,widget_plot_types]) # removed widget_plot_step and widget_binned_pixel_scale_mode
widgets = column(children=[widget_group_one,widget_group_two,widget_group_three],width=dfs.default_toolbar_width)
inputs = row(children=[widgets,tabs],sizing_mode='scale_height')

l = layout([[widget_header],[inputs]])
show(l)

print('Build completed in {} seconds.'.format(time.time()-time_start))