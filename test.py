import pandas as pd
from bokeh.layouts import widgetbox, layout, row
from bokeh.models import ColumnDataSource, CustomJS
from bokeh.plotting import figure, output_file, show
from bokeh.models.widgets import Button, Paragraph

'''
	test thing
		Pandas dataframe --> Pandas ColumnDataSource \ 
												callback --> Placeholder CDS = Pandas CDS (or whatever math you want)
						Placeholder ColumnDataSource /

'''


output_file("panda_test.html")

# my panda import only differs because I messed with my 14d.txt file a bit
panda = pd.read_csv('core/skybackground/14d.txt', sep=',',names=['wavelength','flux14'],skiprows=1)

# ColumnDataSources are the main data-holders for interactive plotting
panda_source = ColumnDataSource(panda)
source = ColumnDataSource(data=dict(x=[], y=[]))

# this one is straight JS, but this can also be done with CustomJS.from_py_func(my_python_function)
callback = CustomJS(args=dict(source=source,panda_source=panda_source), code="""
	var data = source.data;
	var panda_data = panda_source.data;

	data['x'].push(panda_data['wavelength']);
	data['y'].push(panda_data['flux14']);
	data['x'] = data['x'][0]
	data['y'] = data['y'][0]

	console.log(panda_data['wavelength']);
	console.log(data['x']);

	source.change.emit();
	""")

# user interface
button = Button(label="click me", callback=callback)
words = Paragraph(text="""
	The left plot is the pandas data. The right plot contains an empty placeholder line, which is essentially the user's "session."
	When you click the button, the pandas data is operated on (however one pleases) and is pushed to the user's session.
	If you open the browser console (ctrl-shift-i) you'll see the arrays have printed, and match.
	""",width=800)
plot1 = figure(plot_width=400, plot_height=400,x_range=(0,1), y_range=(0,1000000),title='pandas data')
plot2 = figure(plot_width=400, plot_height=400, x_range=(0,1), y_range=(0,1000000),title='placeholder data')
plot_row = row([plot1,plot2])
plot1.line('wavelength','flux14',source=panda_source, line_width=1,line_alpha=0.6, color="blue")
plot2.line('x', 'y', source=source, line_width=1, line_alpha=0.6, color="red")
layout = layout([words,button,plot_row])

show(layout)