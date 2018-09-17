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

base_wavelength = np.arange(start=dfs.default_start_wavelength,stop=dfs.default_stop_wavelength,step=dfs.default_plot_step_step)
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
widget_binned_pixel_scale = RadioButtonGroup(name=dfs.string_binned_pixel_scale_header, labels=dfs.string_binned_pixel_scale_labels,active=0)
# dropdown menus
widget_types_types = Dropdown(label=dfs.string_object_types_types[0],menu=dfs.string_star_types) # TODO: dynamically set these
widget_filters = Dropdown(label=dfs.string_widget_labels[0],menu=dfs.string_filters_menu,width=100)
# text input
widget_mag_input = TextInput(value=dfs.default_magnitude_two,title=dfs.string_suffixes[0].title(),width=75)
widget_redshift = TextInput(value=str(dfs.default_redshift), title=dfs.string_title[3])
#widget_binned_pixel_scale = TextInput(value=dfs.default_binned_pixel_scale,title=dfs.string_binned_pixel_scale_manual) # TODO: hide on default
# sliders
widget_exposure_time = Slider(start=dfs.default_exposure_time_start, end=dfs.default_exposure_time_end,
	value=dfs.default_exposure_time, step=dfs.default_exposure_time_step, title=dfs.string_title[4])
widget_seeing = Slider(start=dfs.default_seeing_start,end=dfs.default_seeing_end,
	value=dfs.default_seeing, step=dfs.default_seeing_step,title=dfs.string_title[5])
widget_slit_width = Slider(start=dfs.default_slit_width_start,end=dfs.default_slit_width_end,
	value=dfs.default_slit_width,step=dfs.default_slit_width_step,title=dfs.string_title[6])
# range sliders
widget_wavelengths = RangeSlider(start=dfs.default_start_wavelength, end=dfs.default_stop_wavelength,
	value=(dfs.default_start_wavelength,dfs.default_stop_wavelength), step=dfs.default_wavelength_step, title=dfs.string_title[8])


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
tabs = Tabs(tabs=[tab0,tab1,tab2,tab3,tab4,tab5,tab6,tab7]) # string_title[10]

# dichroic throughput
cds_dichroic_red = ColumnDataSource(dict(xr=dichro_x,yr=dichro_y1))
cds_dichroic_blue = ColumnDataSource(dict(xb=dichro_x,yb=dichro_y2))
gly_dichroic_red = glyphs.Line(x="xr",y="yr",line_width=dfs.default_line_width,line_color=dfs.default_line_color_red)
gly_dichroic_blue = glyphs.Line(x="xb",y="yb",line_width=dfs.default_line_width,line_color=dfs.default_line_color_blue)
p4.add_glyph(cds_dichroic_red,gly_dichroic_red)
p4.add_glyph(cds_dichroic_blue,gly_dichroic_blue)

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


''' JavaScript functions '''
active_skyfile = []
cb_skyfiles = CustomJS(args=dict(selected=widget_moon_days,activate=active_skyfile), code="""
	activate = selected.active
	""")

resampled_flux = []
cb_res_spect = CustomJS(args=dict(widget=widget_moon_days,skyfiles=coalesced_sky_files[widget_moon_days.active].to_json(),old_spec=base_wavelength,resampled_flux=resampled_flux),code="""	
	var new_spec = skyfiles[widget.active][0] * 1000;
	var spec_flux = skyfiles[widget.active][1] / 1000;
	var ResSpect,
	  splice = [].splice;

	Array.prototype.clone = function() {
	  return this.slice(0);
	  };

	ResSpect = function(new_spec, old_spec, spec_flux) {
	  var end_factor, filter_lhs, filter_width, ref, ref1, ref2, ref3, to_reduce, two_reduce, three_reduce, resampled_flux_err, spec_lhs, spec_width, start, start_factor, stop;
	  spec_lhs = [];
	  spec_width = [];
	  filter_lhs = [];
	  filter_width = [];
	  spec_lhs[0] = old_spec[0] - (old_spec[1] - old_spec[0]) / 2;
	  spec_width[-1] = old_spec[-1] - old_spec[-2];
	  splice.apply(spec_lhs, [1, 0].concat(ref = (old_spec.slice(1, 1) + old_spec.slice(0)) / 2)), ref;
	  splice.apply(spec_width, [0, 0].concat(ref1 = spec_lhs.slice(1, 1) - spec_lhs.slice(0))), ref1;
	  filter_lhs[0] = new_spec[0] - (new_spec[-1] - new_spec[0]) / 2;
	  filter_widths[-1] = new_spec[-1] - new_spec[-2];
	  filter_lhs[-1] = new_spec[-1] + (new_spec[-1] - new_spec[-2]) / 2;
	  splice.apply(filter_lhs, [1, -1].concat(ref2 = (new_spec.slice(1, 1) + new_spec.slice(0)) / 2)), ref2;
	  splice.apply(filter_width, [0, 0].concat(ref3 = filter_lhs.slice(1) - filer_lhs.slice(0, -1))), ref3;
	  if (filter_lhs[0] < spec_lhs[0] || filter_lhs[-1] > spec_lhs[-1]) {
	    console.log("ResSpect: The new wavelengths specified must fall within the range of the old wavelength values.");
	    }
	  resampled_flux = [];
	  if (spec_err === !None) {
	    if (spec_err.shape(!spec_flux.shape)) {
	      console.log("ResSpect: If specified, spec_err must be the same shape as spec_flux.");
	    } else {
	      resampled_flux_err = [];
	      }
	  }
	  start = 0;
	  stop = 0;
	  var j;
	  for (j = 0; j < new_spec.length; j++) {
	    while (spec_lhs[start+1] <= filter_lhs[j]) {
	      start++;
	      }
	    while (spec_lhs[stop+1] <= filter_lhs[j+1]) {
	      stop++;
	      }
	    if (stop === start) {
	      resampled_flux[j] = spec_flux[j] && (spec_errs === !None ? resampled_flux_err[j] = spec_errs[start] : void 0);
	    } else {}        
	        var thing1 = spec_width.slice(start, +(stop + 1) + 1 || 9e9);
	        start_factor = (spec_lhs[start + 1] - filter_lhs[j]) / (spec_lhs[start + 1] - spec_lhs[start]) && (end_factor = (filter_lhs[j + 1] - spec_lhs[stop]) / (spec_lhs[stop + 1] - spec_lhs[stop]) && (spec_widths[start] *= start_factor && (spec_widths[stop] *= end_factor && (resampled_flux[j] = (spec_widths.slice(start, (stop + 1)) * spec_flux.slice(start, (stop + 1))) && (samp_sum += resampled_flux[j] && (!(isNaN(spec_err)) ? void 0 : void 0))))));
	        to_reduce = Math.pow(thing1 * spec_err.slice(start, +(stop + 1) + 1 || 9e9), 2);
	        two_reduce = (spec_width[start..(stop+1)]);
	        three_reduce = Math.pow(thing1 * spec_flux.slice(start, +(stop + 1) + 1 || 9e9), 2);
	        resampled_flux[j] = Math.sqrt(to_reduce.reduceRight((accumulator, currentValue) => accumulator.concat(currentValue)))/three_reduce.reduce((accumulator, currentValue) => accumulator.concat(currentValue));
	        
	        if (spec_err == !None) {
	          resampled_flux_err[j] = Math.sqrt(to_reduce.reduceRight((accumulator, currentValue) => accumulator.concat(currentValue)))/two_reduce.reduce((accumulator, currentValue) => accumulator.concat(currentValue));
	          }
	        spec_width[start] /= start_factor;
	        spec_width[stop] /= end_factor;
	        }
	      }
	  if (spec_err == !None) {
	    return resampled_flux, resampled_flux_err;
	  } else {
	    return resampled_flux;
	    }
	  };
	""")

cb_wavelength = CustomJS(args=dict(p0=p0,p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7,
	slider_wavelengths=widget_wavelengths,slider_plotstep=widget_plot_step,wavelength=wavelength),code="""
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

	p0.x_range.start = a;
	p0.x_range.end = o;
	p0.x_range.change.emit();

	p1.x_range.start = a;
	p1.x_range.end = o;
	p1.x_range.change.emit();

	p2.x_range.start = a;
	p2.x_range.end = o;
	p2.x_range.change.emit();

	p3.x_range.start = a;
	p3.x_range.end = o;
	p3.x_range.change.emit();

	p4.x_range.start = a;
	p4.x_range.end = o;
	p4.x_range.change.emit();

	p5.x_range.start = a;
	p5.x_range.end = o;
	p5.x_range.change.emit();

	p6.x_range.start = a;
	p6.x_range.end = o;
	p6.x_range.change.emit();

	p7.x_range.start = a;
	p7.x_range.end = o;
	p7.x_range.change.emit();
	""")

# observed sky background
cb_sky = CustomJS(args=dict(grating_opt=widget_grating_types.active,skyfile=coalesced_sky_files[widget_moon_days.active],seeing=widget_seeing.value,slit_size=widget_slit_width.value),code="""
	if (grating_opt(0)) {
		delta_lambda_base = 3.73;
	} else {
		delta_lambda_base = 1.40;
	}

	var extension = seeing * slit_size;
	var sky_x = skyfile[0] * 10000;
	var sky_y_ = skyfile[1] / 10000;
	var old_res = sky_x[1] - sky_x[0];
	
	var delta_lambda = delta_lambda_base * slit_size / 0.7;
	var	sigma = delta_lambda / 2.35482;

	""")

# linkages
widget_wavelengths.callback = cb_wavelength
widget_moon_days.callback = cb_res_spect


# final panel building
widget_group_one = widgetbox(children=[widget_telescope_sizes,widget_object_types,widget_galaxy_type])
widget_group_two = layout([[widget_mag_input],[widget_filters,widget_mag_type]])
widget_group_three = widgetbox(children=[widget_grating_types,widget_redshift,widget_exposure_time,widget_seeing,widget_slit_width,widget_moon_days,widget_wavelengths,widget_binned_pixel_scale_mode,widget_binned_pixel_scale,widget_plot_types,widget_plot_step])
widgets = column(children=[widget_group_one,widget_group_two,widget_group_three],width=dfs.default_toolbar_width)
inputs = row(children=[widgets,tabs],sizing_mode='scale_height')
l = layout([[widget_header],[inputs]])

# export
show(l)