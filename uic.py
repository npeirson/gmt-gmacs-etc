#!/usr/bin/python3
'''
	GMACS ETC: UIC
		(user-interface creator)
			builds the etc webpage

'''
import os
import time
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import interpolate # astropy doesn't seem to have pchip integration?
from astropy import units as u
from util import config as cfg
from util import defaults as dfs
from bokeh.plotting import figure
from bokeh.io import output_file,show
from bokeh.models.callbacks import CustomJS
from bokeh.models import ColumnDataSource, glyphs, ranges
from bokeh.layouts import widgetbox, column, row, layout
from bokeh.models.widgets import Dropdown, RadioButtonGroup, CheckboxButtonGroup, Slider, TextInput, RangeSlider, Paragraph, Panel, Tabs, Div

output_file(dfs.string_pagename, title=dfs.site_title, mode='cdn')


'''
	ODC Data Handling
		this part imports data
			performs necessary operations
				and then dies alone like us all

'''
print('GMACS ETC: ODC')
time_start = time.time()

# importing
def dparse(_path,_files):
	_coalesced_data = {}
	for i,_file in tqdm(enumerate(_files),desc='Loading data',ncols=0):
		_full = os.path.join(_path,_file)
		_value = pd.read_csv(_full,delimiter=',')
		_coalesced_data[i] = _value
	return _coalesced_data

coalesced_galaxy_types = dparse(cfg.galaxy_types_path,cfg.galaxy_types_files)
coalesced_star_types = dparse(cfg.star_types_path,cfg.star_types_files)
coalesced_filters = dparse(cfg.filter_path,cfg.filter_files)
coalesced_sky_files = dparse(cfg.skyfiles_path,cfg.skyfiles_files)
coalesced_grating_files = dparse(cfg.efficiency_path,cfg.efficiency_grating_files)
coalesced_ccd_files = dparse(cfg.efficiency_path,cfg.efficiency_ccd_files)
coalesced_dichroic_files = dparse(cfg.dichroic_path,cfg.dichroic_files)
coalesced_atmo_ext_files = dparse(cfg.atmo_ext_path,cfg.atmo_ext_files)

wavelength = np.arange(start=dfs.default_start_wavelength,stop=dfs.default_stop_wavelength,step=dfs.default_plot_step_step)
# efficiency dots
grating_red_1 = np.dot(coalesced_grating_files[0]['a'],10)
grating_red_2 = coalesced_grating_files[0]['b']
#grating_red = interpolate.PchipInterpolator(grating_red_1,grating_red_2,wavelength)
grating_blue_1 = np.dot(coalesced_grating_files[1]['a'],10)
grating_blue_2 = coalesced_grating_files[1]['b']
#grating_blue = interpolate.PchipInterpolator(grating_blue_1,grating_blue_2,wavelength)

ccd_efficiency_red_1 = np.dot(coalesced_ccd_files[0]['a'],10)
ccd_efficiency_red_2 = coalesced_ccd_files[0]['b']
ccd_efficiency_blue_1 = np.dot(coalesced_ccd_files[1]['a'],10)
ccd_efficiency_blue_2 = coalesced_ccd_files[1]['d'] # data came in like this, idk

dichro_x = np.dot(coalesced_dichroic_files[0]['a'],10) # wavelength in Angstroms
dichro_y1 = coalesced_dichroic_files[0]['b'] # reflectivity, blue channel
dichro_y2 = coalesced_dichroic_files[0]['c']

atmo_ext_x = coalesced_atmo_ext_files[0]['a']
atmo_ext_y = coalesced_atmo_ext_files[0]['b']

# divide (n,x,y) into {(n,x),(n,y)}
coalesced_object_x = {}
coalesced_object_y = {}
for i,coalesced_galaxy_type in tqdm(enumerate(coalesced_galaxy_types),desc='Shaping data',ncols=0):
	coalesced_object_x[i] = coalesced_galaxy_types[i].a
	coalesced_object_y[i] = coalesced_galaxy_types[i].b


'''
	ODC User Interface
		manages all interactivity

'''
# radio button groups
widget_telescope_sizes = RadioButtonGroup(labels=dfs.string_telescope_sizes, active=0)
widget_object_types = RadioButtonGroup(labels=dfs.string_object_types, active=dfs.default_object_types)
widget_galaxy_type = Dropdown(label=dfs.string_title[12],menu=dfs.string_galaxy_types)
widget_mag_type = RadioButtonGroup(labels=dfs.string_magnitude_three,active=0)
widget_grating_types = RadioButtonGroup(labels=dfs.string_grating_types, active=0)
widget_moon_days = RadioButtonGroup(labels=dfs.string_moon_days, active=0) # TODO: need title for this, would be string_title[7]
widget_binned_pixel_scale_mode = RadioButtonGroup(labels=dfs.string_binned_pixel_scale_modes, active=0)
# dropdown menus
widget_types_types = Dropdown(label=dfs.string_object_types_types[0],menu=dfs.string_star_types) # TODO: dynamically set these
widget_filters = Dropdown(label=dfs.string_widget_labels[0],menu=dfs.string_filters_menu,width=100)
# text input
widget_mag_input = TextInput(value=dfs.default_magnitude_two,title=dfs.string_suffixes[0].title(),width=75)
widget_redshift = TextInput(value=str(dfs.default_redshift), title=dfs.string_title[3])
widget_seeing = TextInput(value=dfs.default_seeing, title=dfs.string_title[5])
widget_slit_width = TextInput(value=dfs.default_slit_width, title=dfs.string_title[6])
#widget_binned_pixel_scale = TextInput(value=dfs.default_binned_pixel_scale,title=dfs.string_binned_pixel_scale_manual) # TODO: hide on default
widget_binned_pixel_scale = RadioButtonGroup(name=dfs.string_binned_pixel_scale_header, labels=dfs.string_binned_pixel_scale_labels,active=0)

# sliders
widget_exposure_time = Slider(start=dfs.default_exposure_time_start, end=dfs.default_exposure_time_end,
	value=dfs.default_exposure_time, step=1, title=dfs.string_title[4])
widget_wavelengths = RangeSlider(start=dfs.default_start_wavelength, end=dfs.default_stop_wavelength,
	value=(dfs.default_start_wavelength,dfs.default_stop_wavelength), step=1, title=dfs.string_title[8])
# other widgets
widget_header = Div(text='<h1>'+dfs.string_header_one+'</h1><h3>'+dfs.string_header_two+'</h3>',width=500,height=70)
widget_plot_types = CheckboxButtonGroup(labels=dfs.string_plot_types, active=dfs.default_plot_types) # TODO: make double-off == double-on
widget_message_box = Paragraph(text='hello')
widget_plot_step = Slider(start=dfs.default_plot_step[1], end=dfs.default_plot_step[2], value=dfs.default_plot_step[0], step=dfs.default_plot_step_step, title=dfs.string_plot_step)

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

# assemble graph tabs
tab0 = Panel(child=p0,title=dfs.string_calculate_types[0])
tab1 = Panel(child=p1,title=dfs.string_calculate_types[1])
tab2 = Panel(child=p2,title=dfs.string_calculate_types[2])
tab3 = Panel(child=p3,title=dfs.string_calculate_types[3].title())
tab4 = Panel(child=p4,title=dfs.string_calculate_types[4].title())
tab5 = Panel(child=p5,title=dfs.string_calculate_types[5].title())
tab6 = Panel(child=p6,title=dfs.string_calculate_types[6])
tab7 = Panel(child=p7,title=dfs.string_calculate_types[7].title())
tabs = Tabs(tabs=[tab0,tab1,tab2,tab3,tab4,tab5,tab6,tab7]) # string_title[10]


# ccd efficiency
cds_ccd_red = ColumnDataSource(dict(xr=ccd_efficiency_red_1,yr=ccd_efficiency_red_2))
cds_ccd_blue = ColumnDataSource(dict(xb=ccd_efficiency_blue_1,yb=ccd_efficiency_blue_2))
gly_ccd_red = glyphs.Line(x="xr",y="yr",line_width=dfs.default_line_width)
gly_ccd_blue = glyphs.Line(x="xb",y="yb",line_width=dfs.default_line_width)
p6.add_glyph(cds_ccd_red,gly_ccd_red)
p6.add_glyph(cds_ccd_blue,gly_ccd_blue)

# atmospheric extinction
cds_atmo_ext = ColumnDataSource(dict(x=atmo_ext_x, y=atmo_ext_y))
gly_atmo_ext = glyphs.Line(x="x",y="y",line_width=dfs.default_line_width)
p7.add_glyph(cds_atmo_ext,gly_atmo_ext)


''' maybe callbacks go here '''
cb_wavelength = CustomJS(args=dict(fig=p7,slider_wavelengths=widget_wavelengths,slider_plotstep=widget_plot_step,wavelength=wavelength),code="""
	const [a,o] = slider_wavelengths.value;
	var [start,edge] = [a,o];
	var step = slider_plotstep.value;
	var length = edge - start;
	if (length == 1) {
	edge = start;
	start = 0;
	}

	// Validate the edge and step numbers.
	edge = edge || 0;
	step = step || 1;

	// Create the array of numbers, stopping befor the edge.
	for (var ret = []; (edge - start) * step > 0; start += step) {
	ret.push(start);
	}
	wavelength = ret;
	fig.x_range.start = a;
	fig.x_range.end = o;
	console.log(fig.x_range.start);
	console.log(fig.x_range.end);
	fig.x_range.change.emit();
	""")

widget_wavelengths.




''' non-callback plots (no interactivity) '''

# dichroic throughput

# grating throughput



# final panel building
widget_group_one = widgetbox(children=[widget_telescope_sizes,widget_object_types,widget_galaxy_type])
widget_group_two = layout([[widget_mag_input],[widget_filters,widget_mag_type]])
widget_group_three = widgetbox(children=[widget_grating_types,widget_redshift,widget_exposure_time,widget_seeing,widget_slit_width,widget_moon_days,widget_wavelengths,widget_binned_pixel_scale_mode,widget_binned_pixel_scale,widget_plot_types,widget_plot_step])
widgets = column(children=[widget_group_one,widget_group_two,widget_group_three],width=dfs.default_toolbar_width)
inputs = row(children=[widgets,tabs],sizing_mode='scale_height')
l = layout([[widget_header],[inputs]])

# export
show(l)