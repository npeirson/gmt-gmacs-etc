#!/usr/bin/python3

import paths as etpaths

# title, x, y
plot_labels = [('Signal-to-Noise Ratio','Wavelength [\u212b]','SNR px\u207b\u00b9'),
				('Observed Spectrum','Wavelength [\u212b]','Counts px\u207b\u00b9'),
				('Observed Sky Background','Wavelength [\u212b]','Counts px\u207b\u00b9'),
				('Dichroic Throughput','Wavelength [\u212b]','Throughput'),
				('Grating Throughput','Wavelength [\u212b]','Throughput'),
				('CCD Quantum Efficiency','Wavelength [\u212b]','Quantum Efficiency'),
				('Atmospheric Extinction','Wavelength [\u212b]','Extinction')]

widget_headers = ["Telescope Mode", # 0
					"Object Type", # 1
					"Stellar Classification", # 2
					"Galactic Classifications", # 3
					"Magnitude System", # 4
					"Magnitude", # 5
					"Filter", # 6
					"Grating", # 7
					"Days since/until new moon", # 8
					"Pixel Binning", # 9
					"Redshift", # 10
					"Seeing", # 11
					"Slit Width", # 12
					"Exposure Time", # 13
					"Spectral Range", # 14
					"Include Noise", # 15
					"Active Channels"] # 16

widget_names = ['widget_telescope','widget_object_type','widget_star_type',
				'widget_galaxy_type','widget_mag_sys','widget_mag',
				'widget_filter','widget_grating','widget_moon',
				'widget_binning','widget_redshift','widget_seeing',
				'widget_slit','widget_time','widget_wavelength',
				'widget_withnoise','widget_channels'] # matched w/ header nums

telescope_sizes = ["First light","Full Size"]
object_types = ["Stellar","Galactic"]

star_types = [name[:-4] for name in etpaths.stellar_files]
galaxy_types = etpaths.galaxy_files # these files lack suffixes
star_types_tup = [(i,i) for i in star_types] # unless there's a need to distinguish...
galaxy_types_tup = [(i,i) for i in galaxy_types]
filters_tup = [(name[:-4],etpaths.filter_files[i]) for i,name in enumerate(etpaths.filter_files)]


mag_sys_opts = ['Vega','AB']
grating_opts = ['Low Resolution','High Resolution']
filter_opts = [name[:-4] for name in etpaths.filter_files]
moon_opts = ['0','3','7','10','14']
bin_opts = ['1x1','2x2','3x3','4x4']
noise_opts = ['Off','On']
channels = ['Blue','Red']

