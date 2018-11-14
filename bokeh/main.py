from os.path import dirname, join

import numpy as np

from bokeh.plotting import figure
from bokeh.layouts import layout,widgetbox,column,row
from bokeh.models import ColumnDataSource,glyphs
from bokeh.io import curdoc
from bokeh.models.widgets import Dropdown,RadioButtonGroup,CheckboxButtonGroup,Slider,RangeSlider,Tabs,Panel,Div
from bokeh.embed import components # for loading stuff, todo

import slim as etslim # reduced python package for gui versions
import strings as stc
import defaults as dfs


# create widgets
widget_telescope_size = RadioButtonGroup(labels=stc.telescope_sizes,active=0,name=stc.widget_names[0])
widget_object_type = RadioButtonGroup(labels=stc.object_types,active=0,name=stc.widget_names[1])
widget_star_type = Dropdown(default_value=stc.star_types_tup[4][0],label=stc.widget_headers[2],menu=stc.star_types_tup,name=stc.widget_names[2])
widget_galaxy_type = Dropdown(default_value=stc.galaxy_types_tup[0][0],label=stc.widget_headers[3],menu=stc.galaxy_types_tup,name=stc.widget_names[3])
widget_mag_sys = RadioButtonGroup(labels=stc.mag_sys_opts,active=1,name=stc.widget_names[4])
widget_mag = Slider(start=(-27),end=(32),value=(25),step=(0.1),title=stc.widget_headers[5],name=stc.widget_names[5])
widget_filter = Dropdown(default_value=stc.filters_tup[7][0],label=stc.widget_headers[6],menu=stc.filters_tup,name=stc.widget_names[6])
widget_grating = RadioButtonGroup(labels=stc.grating_opts,active=0,name=stc.widget_names[7])
widget_moon = RadioButtonGroup(labels=stc.moon_opts,active=0,name=stc.widget_names[8])
widget_binning = RadioButtonGroup(labels=stc.bin_opts,active=1,name=stc.widget_names[9])
widget_redshift = Slider(start=(0),end=(13.37),value=(0),step=(0.01),title=stc.widget_headers[10],name=stc.widget_names[10])
widget_seeing = Slider(start=(0),end=(20),value=(0.5),step=(0.1),title=stc.widget_headers[11],name=stc.widget_names[11])
widget_slit = Slider(start=(0),end=(10),value=(0.5),step=(0.05),title=stc.widget_headers[12],name=stc.widget_names[12])
widget_time = Slider(start=(0),end=(21600),value=(2600),step=(10),title=stc.widget_headers[13],name=stc.widget_names[13])
widget_wavelength = RangeSlider(start=dfs.wavelength_limits[0], end=dfs.wavelength_limits[1],
    value=((dfs.wavelength_limits[0]+1400),(dfs.wavelength_limits[1]-1400)),step=(10),title=stc.widget_headers[14],name=stc.widget_names[14])
widget_withnoise = RadioButtonGroup(labels=stc.noise_opts,active=1,name=stc.widget_names[15])
widget_channels = CheckboxButtonGroup(labels=stc.channels, active=[0,1], name=stc.widget_names[16])

widget_header = Div(text='<h1>'+stc.header1+'</h1><h3>'+stc.header2+'</h3>',width=500,height=70)

# create figures
p0 = figure(plot_width=dfs.plot_dims[0], plot_height=dfs.plot_dims[1],sizing_mode=dfs.plot_sizing_mode,
            x_axis_label=stc.plot_labels[0][0],y_axis_label=stc.plot_labels[0][1])
p1 = figure(plot_width=dfs.plot_dims[0], plot_height=dfs.plot_dims[1],sizing_mode=dfs.plot_sizing_mode,
            x_axis_label=stc.plot_labels[1][0],y_axis_label=stc.plot_labels[1][1])
p2 = figure(plot_width=dfs.plot_dims[0], plot_height=dfs.plot_dims[1],sizing_mode=dfs.plot_sizing_mode,
            x_axis_label=stc.plot_labels[2][0],y_axis_label=stc.plot_labels[2][1])
p3 = figure(plot_width=dfs.plot_dims[0], plot_height=dfs.plot_dims[1],sizing_mode=dfs.plot_sizing_mode,
            x_axis_label=stc.plot_labels[3][0],y_axis_label=stc.plot_labels[3][1])
p4 = figure(plot_width=dfs.plot_dims[0], plot_height=dfs.plot_dims[1],sizing_mode=dfs.plot_sizing_mode,
            x_axis_label=stc.plot_labels[4][0],y_axis_label=stc.plot_labels[4][1])
p5 = figure(plot_width=dfs.plot_dims[0], plot_height=dfs.plot_dims[1],sizing_mode=dfs.plot_sizing_mode,
            x_axis_label=stc.plot_labels[5][0],y_axis_label=stc.plot_labels[5][1])
p6 = figure(plot_width=dfs.plot_dims[0], plot_height=dfs.plot_dims[1],sizing_mode=dfs.plot_sizing_mode,
            x_axis_label=stc.plot_labels[6][0],y_axis_label=stc.plot_labels[6][1])
figures = [p0,p1,p2,p3,p4,p5,p6] # group figures

# create figure tabs
tab0 = Panel(child=p0,title=stc.plot_labels[0][0])
tab1 = Panel(child=p1,title=stc.plot_labels[1][0])
tab2 = Panel(child=p2,title=stc.plot_labels[2][0])
tab3 = Panel(child=p3,title=stc.plot_labels[3][0])
tab4 = Panel(child=p4,title=stc.plot_labels[4][0])
tab5 = Panel(child=p5,title=stc.plot_labels[5][0])
tab6 = Panel(child=p6,title=stc.plot_labels[6][0])
widget_tabs = Tabs(tabs=[tab0,tab1,tab2,tab3,tab4,tab5,tab6],name='widget_tabs')

# group widgets for initialization (not layout)
widgets_with_active = [widget_telescope_size,widget_object_type,widget_mag_sys,
                        widget_grating,widget_moon,widget_binning,widget_withnoise,widget_channels]
widgets_with_values = [widget_star_type,widget_galaxy_type,widget_filter,widget_mag,
                        widget_redshift,widget_seeing,widget_slit,widget_time]
widgets_coalesced = np.append(np.append(np.append(widgets_with_active,widgets_with_values),widget_wavelength),widget_tabs)

# create glyphs and their sources
cds_blue = ColumnDataSource(data=dict(xb=[],yb=[]))
cds_red = ColumnDataSource(data=dict(xr=[], yr=[]))
gly_blue = glyphs.Line(x='xb',y='yb',line_color='blue')
gly_red = glyphs.Line(x='xr',y='yr',line_color='red')
for p in figures:
    p.add_glyph(cds_blue,gly_blue)
    p.add_glyph(cds_red,gly_red)
val_pack = dict()
names = [widgets_coalesced[i].name[7:] for i in range(len(widgets_coalesced))]
[val_pack.update({names[i]:widgets_coalesced[i]}) for i in range(widgets_coalesced.shape[0])] # comprehend this, yo
sess = etslim.session(val_pack) # create an etc session object with initial values

def update_bkh(caller):
    plot_x,plot_yb,plot_yr = sess.update(caller)
    if 0 in widget_channels.active:
        cds_blue.data['xb'] = plot_x
        cds_blue.data['yb'] = plot_yb
    if 1 in widget_channels.active:
        cds_red.data['xr'] = plot_x
        cds_red.data['yr'] = plot_yr
    else: # crashless catch-all
        cds_blue.data['xb'] = plot_x
        cds_blue.data['yb'] = plot_yb
        cds_red.data['xr'] = plot_x
        cds_red.data['yr'] = plot_yr

# link callbacks
for i,widge in enumerate(widgets_with_values):
    print(names[i])
    widge.on_change('value', lambda attr, old, new: update_bkh(names[i]))
for i,widge in enumerate(widgets_with_active):
    print(names[i])
    widge.on_change('active', lambda attr, old, new: update_bkh(names[i]))

sizing_mode = 'scale_both'

widget_group_one = widgetbox(children=[widget_telescope_size,widget_object_type,widget_star_type,widget_galaxy_type])
widget_group_two = layout([[widget_mag],[widget_filter,widget_mag_sys]])
widget_group_three = widgetbox(children=[widget_grating,widget_redshift,widget_time,widget_seeing,widget_slit,widget_moon,widget_wavelength,widget_binning,widget_channels])
widgets = column(children=[widget_group_one,widget_group_two,widget_group_three],width=dfs.toolbar_width)
inputs = row(children=[widgets,widget_tabs],sizing_mode='scale_height')

l = layout([[widget_header],[inputs]])

#update()  # initial load of the data

curdoc().add_root(l)
curdoc().title = "GMACS ETC"