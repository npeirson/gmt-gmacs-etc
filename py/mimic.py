import math
import numpy as np
import values as edl
from spectres import spectres
from astropy import constants as const
from astropy.convolution import convolve_fft
from astropy.stats import gaussian_sigma_to_fwhm

string_prefix = '[ etc ] :'
coating_eff = 0.8*0.98**14

stellar_type_keys = ['star','stars','stellar',0]
galactic_type_keys = ['galaxies','galaxy','galactic',1]
grating_opt_keys = ['low',0,1.4,'high',1,3.73]
moon_days_keys = [0,3,7,10,14]
telescope_mode_keys = [0,4,'first','first light',1,7,'full','full size']

if isinstance(object_type,str):
	object_type = object_type.lower()
if isinstance(grating_opt,str):
	grating_opt = grating_opt.lower()
if isinstance(telescope_mode,str):
	telescope_mode = telescope_mode.lower()

if object_type in stellar_type_keys:
	type_keys = edl.stellar_files + [name[:-4] for name in in edl.stellar_files]
elif object_type in galactic_type_keys:
	type_keys = edl.galaxy_files + [name[:-4] for name in in edl.stellar_files]
else:
	raise ValueError("{} Invalid object type ({})".format(string_prefix,object_type))

if grating_opt in grating_opt_keys[:2]: # low resolution
	delta_lambda_default = edl.dld[1]
elif grating_opt in grating_opt_keys[3:]:
	delta_lambda_default = edl.dld[0]
else:
	raise ValueError("{} Invalid grating option ({})".format(string_prefix,grating_opt))

delta_lambda = delta_lambda_default * slit_size / 0.7

if moon_days in moon_days_keys:
	skyfile = str(moon_days) + 'd.txt'
else:
	raise ValueError('{} Invalid number of days since new moon ({})'.format(string_prefix,moon_days))

if telescope_mode in telescope_mode_keys[:3]:
	area = edl.area[0]
elif telescope_mode in telescope_mode_keys[4:]:
	area = edl.area[1]
else:
	raise ValueError('{} Invalid telescope mode ({})'.format(string_prefix,telescope_mode))

if ((wavelength[0] > 3200) and (wavelength[-1] < 12000)) or ((wavelength[0] > 320) and (wavelength[-1] < 1200)):
	plot_step = wavelength[2] - wavelength[1]
elif wavelength is 'default':
	plot_step = delta_lambda / 3
else:
	raise ValueError('{} Invalid wavelength range ({}-{})'.format(string_prefix,wavelength[0],wavelength[-1]))

def mag_cal(wavelength,filter_opt,mag_sys_opt,object_type,redshift,mag):
	'''
		sample inputs:
		wavelength = np.ndarray
		filter_opt = 'r'
		mag_sys_opt = 'vega'
		object_type = 'a0v'
		redshift = 0
		mag = 25
	'''
	lambda_min = min(selected_filter[0])
	lambda_max = max(selected_filter[0])
	lambda_A = np.arange(lambda_min,lambda_max,plot_step)

	object_x = np.dot(selected_object[0],(1+redshift))
	object_y = selected_object[1]
	flux_A = spectres(lambda_A,object_x,object_y)

	trans = spectres(lambda_A,filter_x,filter_y)

	extinction = spectres(lambda_A,atmo_ext[0],atmo_ext[1])

	flux = np.dot(flux_A,1e10)
	_lambda = np.dot(lambda_A,1e10)

	num_zeros = 0
	for lux in flux:
		if (lux is None) or (lux is 0):
			num_zeros += 1
			lux = 0

	if (num_zeros >= (flux.shape[0]/5)):
		if (num_zeros == flux.shape[0]):
			print('No flux in this bandpass!')
			output_flux = np.zeros(wavelength.shape[0])
		else:
			percent_zeros = (num_zeros / flux.shape[0]) * 100
			print('{}% of this bandpass has zero flux'.format(percent_zeros))

	if (mag_sys_opt == 'vega'):
		flux_vega = spectres(wavelength,vega[0],vega[1])
		mag_model = -2.5 * np.log10(np.divide(fsum(flux.dot(trans).dot(extinction).dot(_lambda)),fsum(flux_vega.dot(trans).dot(extinction).dot(_lambda)))) + 0.03
	elif (mag_sys_opt == 'ab'):
		mag_model = -48.60 - 2.5 * np.log10(np.divide(fsum(flux.dot(trans).dot(extinction).dot(_lambda)),fsum(trans.dot(_lambda).dot(extinction).dot(np.divide(const.c,(np.square(_lambda)))))))
	else:
		print('Invalid mag_sys_opt!')

	del_mag = mag_input - mag_model
	output_lambda = object_x
	output_flux = np.dot(object_y,np.negative(np.divide(del_mag/2.5)))
	return output_lambda, output_flux

star_x, star_y = mag_cal(wavelength,filter_opt,mag_sys_opt,object_type,redshift,mag)

old_res = star_x[2] - star_x[1]
if (old_res < plot_step):
	_sigma = delta_lambda / gaussian_sigma_to_fwhm
	_x = np.arange((-5*_sigma),(5*_sigma),old_res) # lambda function this?
	funx = lambda x: (1/(_sigma*np.sqrt(2*math.pi)))*np.exp(np.divide(np.negative(np.square(x)),(np.dot(np.square(_sigma),2))))
	degrade = funx(_x)/np.trapz(funx(_x))
	star_y = convolve_fft(star_y,degrade)
else:
	flux_y = spectres(star_x,star_y,wavelength)

flux = np.dot(flux_y,1e-03)
power = np.dot(np.dot(flux,area),np.dot(exp_time,plot_step))
counts = np.divide(np.divide(power,np.divide((const.h.value * const.c.value),wavelength)),1e10)


''' subtract light lost to various components '''

# seeing
_sigma = seeing / gaussian_sigma_to_fwhm
percent_u,percent_err_u = quad(funx,(-self.slit_size/2),(self.slit_size/2))
percent_l,percent_err_l = quad(funx,(-self.seeing/2),(self.seeing/2))
percent = percent_u * percent_l # can use error if you add it later...
extension = seeing * slit_size

# sky background
sky_x = np.dot(sky_background[0],1e4)
sky_y = np.divide(sky_background[1],1e4)
old_res = sky_x[2] - sky_x[1]
_sigma = delta_lambda / gaussian_sigma_to_fwhm
_x = np.arange((-5*_sigma),(5*_sigma),old_res)
degrade = funx(_x)/np.trapz(funx(_x))
sky_y = convolve_fft(sky_y,degrade)
sky_flux = spectres(wavelength,sky_x,sky_y)
counts_noise = np.dot(np.dot(sky_flux,extension),(area*exp_time*plot_step))

# dichroic
dichro_x = dichro[0] * 10
dichro_y1 = dichro[1]
dichro_y2 = dichro[2]
blue_dichro = spectres(wavelength,dichro_x,dichro_y1)
red_dichro = spectres(wavelength,dichro_x,dichro_y2)

# grating
blue_grating = spectres(wavelength,(grating1[0]*10),grating1[1])
red_grating = spectres(wavelength,(grating2[0]*10),grating2[1])

# ccd
blue_ccd = spectres(wavelength,(ccd1[0]*10),ccd1[1])
red_ccd = spectres(wavelength,(ccd2[0]*10),ccd2[1])

# atmospheric extinction
atmo_ext_x = atmo_ext[0]
atmo_ext_y = atmo_ext[1]

# mirror bouncing, apply twice
mirror = spectres(wavelength,mirror_file[0]*10,mirror_file[1])


''' calculations '''

# signal
blue_total_eff = np.dot(np.dot(blue_dichro, blue_grating), (blue_ccd * (coating_eff * extinction)))
red_total_eff = np.dot(np.dot(red_dichro, red_grating), (red_ccd * (coating_eff * extinction)))
blue_signal = np.dot((counts * percent), blue_total_eff)
red_signal = np.dot((counts * percent), red_total_eff)

# noise
blue_total_eff_noise = np.dot(np.dot(blue_dichro,blue_grating),(blue_ccd * coating_eff))
red_total_eff_noise = np.dot(np.dot(red_dichro,red_grating),(red_ccd * coating_eff))
blue_noise = np.dot(counts_noise * blue_total_eff_noise)
red_noise = np.dot(counts_noise * red_total_eff_noise)

snr_blue = np.divide(blue_signal,np.sqrt(blue_signal + blue_noise + np.square(readnoise)))
snr_red = np.divide(red_signal,np.sqrt(red_signal + red_noise + np.square(readnoise)))

sigma_blue = np.sqrt(blue_signal + blue_noise)
sigma_red = np.sqrt(red_signal + red_noise)

error_blue = np.random.normal(loc=0, scale=(np.sqrt(blue_signal + blue_noise)),size=len(blue_snr))
error_red = np.random.normal(loc=0, scale=(np.sqrt(red_signal + red_noise)),size=len(red_snr))


''' pre-plotting '''

thing_keys = ['snr','obv_spec','sky_background','dichroic_throughput','grating_throughput','ccd_qe','atmospheric_extinction']
plot_x = wavelength
if thing in thing_keys:
	if (thing.lower() is thing_keys[0]):
		title = 'Signal-to-Noise Ratio'
		labels = ['Angstrom','Signal-to-Noise Ratio']
		plot_y_blue = blue_snr
		plot_y_red = red_snr
	elif (thing.lower() is thing_keys[1]):
		labels = ['Angstrom','Counts (#) per pixel']
		if noise:
			title = 'Observed Spectrum (with noise)'
			plot_y_blue = np.add(blue_signal,error_blue)
			plot_y_red = np.add(red_signal,error_red)
		else:
			title = 'Observed Spectrum (without noise)'
			plot_y_blue = blue_signal
			plot_y_red = red_signal
	elif (thing.lower() is thing_keys[2]):
		title = 'Observed Sky Background'
		labels = ['Angstrom','Counts (#) per pixel']
		plot_y_blue = blue_noise
		plot_y_red = red_noise
	elif (thing.lower() is thing_keys[3]):
		title = 'Dichroic Throughput'
		labels = ['Angstrom','Throughput']
		plot_y_blue = blue_dichro
		plot_y_red = red_dichro
	elif (thing.lower() is thing_keys[4]):
		title = 'Grating Throughput'
		labels = ['Angstrom','Throughput']
		plot_y_blue = blue_grating
		plot_y_red = red_grating
	elif (thing.lower() is thing_keys[5]):
		title = 'CCD Quantum Efficiency'
		labels = ['Angstrom','QE']
		plot_y_blue = blue_ccd
		plot_y_red = red_ccd
	elif (thing.lower() is thing_keys[6]):
		title = 'Atmospheric Extinction'
		labels = ['Angstrom','Throughput']
		plot_y_blue = spectres(wavelength,atmo_ext_x,atmo_ext_y)
		plot_y_red = plot_y_blue
	else:
		raise ValueError('{} Invalid thing ({})'.format(string_prefix,thing))
