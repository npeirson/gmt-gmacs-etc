# this uses a caller singleton and enumerated caller-association groups to update only necessary data.
# the else statement catches errors and wavelength changes.

import values as edl

class simulate:

	def __init__(self):
		pass # do init stuff!


	''' tier 0 '''

	def refresh(self,caller,tabs=tabs): # default cb_obj as immutable
		if caller in edl.signal_components:
			recalculate_signal(caller)
		if caller in edl.noise_components:
			recalculate_noise(caller)
		else:
			recalculate_signal(caller)
			recalculate_noise(caller)
		if (tabs.active == 0): # snr
			recalculate_snr(caller)
		elif (tabs.active == 1): # obs spec no noise
			pass # just send signal back :)
		elif (tabs.active == 2): # obs spec noise
			recalculate_error(caller)
		elif (tabs.active == 3): #obs sky background
			pass # send noise
		return the_stuff


	def  recalculate_snr(self,caller):
		pass


	def recalculate_error(self,caller):
		pass


	''' tier 1 '''

	def recalculate_signal(self,caller):
		if caller in edl.counts_components:
			pass # recalc counts stuff
		elif caller in edl.percent_components:
			pass # recalc percent stuff
		elif caller in edl.efficiency_components:
			pass # recalc efficiency stuff
		else:
			pass # recalc all stuff


	def recalculate_noise(self,caller):
		if caller in edl.counts_noise_components:
			pass # recalc counts noise stuff
		elif caller in edl.efficiency_noise_components:
			pass # recalc total efficiency noise stuff


	''' tier 2 '''

	def recalculate_counts(self,caller):
		pass # recalc counts stuff


	def recalculate_counts_noise(self,caller):
		pass # recalc counts noise stuff


	def recalculate_percent(self,caller):
		pass # recalc percent stuff


	def recalculate_efficiency(self,caller):
		pass # recalc efficiency stuff

	def recalculalte_efficiency_noise(self,caller):
		pass # recalc efficienccy noise stuff

