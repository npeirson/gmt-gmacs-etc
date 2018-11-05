import math
import json
import numpy as np
import values as edl
import datahandler as dh
from spectres import spectres
from scipy import interpolate
from astropy import constants as const
from astropy.convolution import convolve, convolve_fft
from astropy.stats import gaussian_sigma_to_fwhm


''' stuff to move '''
galaxyfiles = [[dh.sb1_x,dh.sb1_y],[dh.sb2_x,dh.sb2_y],[dh.sb3_x,dh.sb3_y],[dh.sb4_x,dh.sb4_y],[dh.sb5_x,dh.sb5_y],
				[dh.sb6_x,dh.sb6_y],[dh.s0_x,dh.s0_y],[dh.sa_x,dh.sa_y],[dh.sb_x,dh.sb_y],[dh.sc_x,dh.sc_y],
				[dh.bulge_x,dh.bulge_y],[dh.ellipticals_x,dh.ellipticals_y],[dh.lbg_all_flam_x,dh.lbg_all_flam_y]]
starfiles = [[dh.star_o5v_x,dh.star_o5v_y],[dh.star_b0v_x,dh.star_b0v_y],[dh.star_b57v_x,dh.star_b57v_y],
			[dh.star_a0v_x,dh.star_a0v_y],[dh.star_a5v_x,dh.star_a5v_y],[dh.star_f0v_x,dh.star_f0v_y],
			[dh.star_g0v_x,dh.star_g0v_y],[dh.star_g5v_x,dh.star_g5v_y],[dh.star_k0v_x,dh.star_k5v_y],
			[dh.star_m0v_x,dh.star_m0v_y],[dh.star_m5v_x,dh.star_m5v_y]]
filterfiles = [[dh.filter_photonux_x,dh.filter_photonux_y],[dh.filter_photonb_x,dh.filter_photonb_y],
				[dh.filter_photonv_x,dh.filter_photonv_y],[dh.filter_photonr_x,dh.filter_photonr_y],
				[dh.filter_photoni_x,dh.filter_photoni_y],[dh.filter_u_x,dh.filter_u_y],[dh.filter_g_x,dh.filter_g_y],
				[dh.filter_r_x,dh.filter_r_y],[dh.filter_i_x,dh.filter_i_y],[dh.filter_z_x,dh.filter_z_y]]
skyfiles = [[dh.skyfile_00d_x,dh.skyfile_00d_y],[dh.skyfile_03d_x,dh.skyfile_03d_y],[dh.skyfile_07d_x,dh.skyfile_07d_y],
			[dh.skyfile_10d_x,dh.skyfile_10d_y],[dh.skyfile_14d_x,dh.skyfile_14d_y]]
atmo_ext_x = dh.atmo_ext_x
atmo_ext_y = dh.atmo_ext_y


''' constants '''

object_type = 'a5v'
filter_index = -3
mag_sys_opt = 'ab'
magnitude = 25
redshift = 0
seeing = 0.5
slit_size = 0.5
moon_days = 0
grating_opt = 0
telescope_mode = 'first'
exp_time = 3600

wavelength = np.arange(3200,10360,10)
channel = 'red'

string_prefix = '[ etc ] :'
coating_eff = 0.8*0.98**14

stellar_keys = [filename[:-4] for filename in edl.stellar_files]
galactic_keys = edl.galaxy_files
filter_keys = [filename[:-4] for filename in edl.filter_files]
grating_opt_keys = ['low',0,1.4,'high',1,3.73]
moon_days_keys = [0,3,7,10,14]
telescope_mode_keys = [0,4,'first','first light',1,7,'full','full size']
selected_filter = filterfiles[filter_index]

''' input handlding '''

if isinstance(object_type,str):
	object_type = object_type.lower()
if isinstance(grating_opt,str):
	grating_opt = grating_opt.lower()
if isinstance(telescope_mode,str):
	telescope_mode = telescope_mode.lower()

if object_type in stellar_keys:
	index_of = [i for i,name in enumerate(stellar_keys) if object_type in name][0]
	object_type = starfiles[index_of]
elif object_type in galactic_keys:
	index_of = [i for i,name in enumerate(galactic_keys) if object_type in name][0]
	object_type = galaxyfiles[index_of]
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
	sky_background = skyfiles[(int(np.where(np.asarray(moon_days_keys)==moon_days)[0]))]
else:
	raise ValueError('{} Invalid number of days since new moon ({})'.format(string_prefix,moon_days))

if telescope_mode in telescope_mode_keys[:3]:
	area = edl.area[0]
elif telescope_mode in telescope_mode_keys[4:]:
	area = edl.area[1]
else:
	raise ValueError('{} Invalid telescope mode ({})'.format(string_prefix,telescope_mode))

plot_step = wavelength[2] - wavelength[1]
if ((wavelength[0] > (3200-plot_step)) and (wavelength[-1] < (10360+plot_step))) or ((wavelength[0] > (320-plot_step)) and (wavelength[-1] < (1036+plot_step))): # picks up either unit
	plot_step = wavelength[2] - wavelength[1]
elif wavelength is 'default':
	plot_step = delta_lambda / 3
else:
	raise ValueError('{} Invalid wavelength range ({}-{})'.format(string_prefix,wavelength[0],wavelength[-1]))


''' mag_cal '''

def mag_cal(wavelength,selected_filter,mag_sys_opt,object_type,redshift,mag):
	'''
		sample inputs:
		wavelength = np.ndarray
		selected_filter = 'r'
		mag_sys_opt = 'vega'
		object_type = 'a0v'
		redshift = 0
		mag = 25
	'''
	filter_min = min(selected_filter[0])
	filter_max = max(selected_filter[0])

	if (filter_min > wavelength[0]):
		lambda_min = filter_min
	elif (filter_min == wavelength[0]):
		filter_min = selected_filter[int(np.where(selected_filter[0] > wavelength[0])[0])]
	else:
		lambda_min = wavelength[0]

	if (filter_max < wavelength[-1]):
		lambda_max = filter_max
	elif (filter_max == wavelength[-1]):
		filter_max = selected_filter[int(np.where(selected_filter[0] < wavelength[-1])[-1])]
	else:
		lambda_max = wavelength[-1]

	plot_step = wavelength[2] - wavelength[1]
	lambda_A = np.arange(lambda_min,lambda_max,plot_step)

	object_x = object_type[0] * (1+redshift)
	object_y = object_type[1]
	flux_A = spectres(lambda_A,object_x,object_y)

	lambda_A[0] = lambda_A[0] + plot_step
	lambda_A[-1] = lambda_A[-1] - plot_step

	trans = spectres(lambda_A,selected_filter[0],selected_filter[1])

	extinction = spectres(lambda_A,atmo_ext_x,atmo_ext_y)

	flux = flux_A * 1e10
	_lambda = lambda_A * 1e10

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


	print(wavelength.shape[0]) # 716
	print(extinction.shape[0]) # 185
	print(_lambda.shape[0])	# 185
	print(trans.shape[0]) # 75

	if (mag_sys_opt == 'vega'):
		flux_vega = spectres(wavelength,vega[0],vega[1])
		mag_model = -2.5 * np.log10(np.divide(math.fsum(flux.dot(extinction).dot(_lambda) * trans[1]),math.fsum(flux_vega.dot(trans).dot(_lambda) * extinction[1]))) + 0.03
	elif (mag_sys_opt == 'ab'):
		#mag_model = -48.60 - 2.5 * np.log10(np.divide(math.fsum(np.dot(flux,trans).dot(extinction).dot(_lambda)),math.fsum(trans.dot(_lambda) * extinction[1]).dot(const.c.value/(np.square(_lambda)))))
		mag_model = -48.6 - 2.5 * np.log10((np.dot(np.dot(flux,trans),np.dot(extinction,_lambda))) / (math.fsum(np.dot((np.dot(trans,_lambda) * extinction[1]),(const.c.value/np.square(_lambda))))))
	else:
		print('Invalid mag_sys_opt!')

	del_mag = mag - mag_model
	output_lambda = object_x
	output_flux = np.dot(object_y,np.negative(del_mag/2.5))
	return output_lambda, output_flux


star_x, star_y = mag_cal(wavelength=wavelength,selected_filter=selected_filter,mag_sys_opt=mag_sys_opt,object_type=object_type,redshift=redshift,mag=magnitude)
print(star_x.shape[0])
print(star_y.shape[0])
old_res = star_x[2] - star_x[1]
if (old_res < plot_step):
	flux_y = spectres(wavelength,star_x,star_y)
else:
	flux_y = spectres(wavelength,star_x,star_y)

flux = np.dot(flux_y,1e-03)
power = np.dot(np.dot(flux,area),np.dot(exp_time,plot_step))
counts = np.divide(np.divide(power,np.divide((const.h.value * const.c.value),wavelength)),1e10)


''' subtract light lost to various components '''

# seeing
_sigma = seeing / gaussian_sigma_to_fwhm
percent_u,percent_err_u = np.trapz(funx,(-self.slit_size/2),(self.slit_size/2))
percent_l,percent_err_l = np.trapz(funx,(-self.seeing/2),(self.seeing/2))
percent = percent_u * percent_l # can use error if you add it later...
extension = seeing * slit_size

# sky background
sky_x = sky_background[0] * 1e4
sky_y = sky_background[1] / 1e4
old_res = sky_x[2] - sky_x[1]
_sigma = delta_lambda / gaussian_sigma_to_fwhm
_x = np.arange((-5*_sigma),(5*_sigma),old_res)
degrade = funx(_x)/np.trapz(funx(_x))
sky_y = convolve_fft(sky_y,degrade)
sky_flux = spectres(wavelength,sky_x,sky_y)
counts_noise = np.dot(np.dot(sky_flux,extension),(area*exp_time*plot_step))

# dichroic
#dichro_x = dichro[0] * 10
#dichro_y1 = dichro[1]
#dichro_y2 = dichro[2]
if (channel is 'blue') or (channel is 'both'):
	blue_dichro = spectres(wavelength,dichro_x,dichro_y1)
if (channel is 'red') or (channel is 'both'):
	red_dichro = spectres(wavelength,dichro_x,dichro_y2)

# grating
if (channel is 'blue') or (channel is 'both'):
	blue_grating = spectres(wavelength,(grating1[0]*10),grating1[1])
if (channel is 'red') or (channel is 'both'):
	red_grating = spectres(wavelength,(grating2[0]*10),grating2[1])

# ccd
if (channel is 'blue') or (channel is 'both'):
	blue_ccd = spectres(wavelength,(ccd1[0]*10),ccd1[1])
if (channel is 'red') or (channel is 'both'):
	red_ccd = spectres(wavelength,(ccd2[0]*10),ccd2[1])

# atmospheric extinction
# should already be globals
# atmo_ext_x = atmo_ext[0]
# atmo_ext_y = atmo_ext[1]

# mirror bouncing, apply twice
mirror = spectres(wavelength,mirror_file[0]*10,mirror_file[1])

# read noise
spectral_resolution = math.ceil((slit_size/(0.7/12))/2)*2 #px (ceil()/2)*2 to round up to next even integer
spatial_resolution = math.ceil((seeing/(0.7/12))/2)*2 #px (ceil()/2)*2 to round up to next even integer 
extent = seeing * slit_size
npix = extent/(0.7/12)**2
try:
	isinstance(bin_size,int)
except:
	bin_size = edl.bin_options_int[edl.bin_options_default_index] # default 2x2 binning

rn = edl.rn_default
if (bin_size > 0) and (bin_size < 5):
	print('[ info ] : Pixel binning: ({}x{})'.format(bin_size,bin_size))
	readnoise = math.ceil(self.rn * spectral_resolution * spatial_resolution / (bin_size**2))
	print('[ info ] : Extent: {} arcsec^2\n[ info ] : num pixels: {} px\n[ info ] : spectral resolution: {} px\n[ info ] : spatial resolution: {} px'.format(self.extent,int(math.ceil(npix)),spectral_resolution,spatial_resolution))
else:
	raise ValueError('{} Invalid pixel binning option ({})'.format(string_prefix,bin_size))


''' calculations '''

# signal
if (channel is 'blue') or (channel is 'both'):
	blue_total_eff = np.dot(np.dot(blue_dichro,blue_grating),np.dot((blue_ccd * (coating_eff * extinction)),np.square(mirror)))
	blue_signal = np.dot((counts * percent), blue_total_eff)
if (channel is 'red') or (channel is 'both'):
	red_total_eff = np.dot(np.dot(red_dichro,red_grating),np.dot((red_ccd * (coating_eff * extinction)),np.square(mirror)))
	red_signal = np.dot((counts * percent), red_total_eff)

# noise
if (channel is 'blue') or (channel is 'both'):
	blue_total_eff_noise = np.dot((np.dot(np.dot(blue_dichro,blue_grating)),(np.dot(blue_ccd,np.square(mirror)) * coating_eff)))
	blue_noise = np.dot(counts_noise * blue_total_eff_noise)
if (channel is 'red') or (channel is 'both'):
	red_total_eff_noise = np.dot((np.dot(np.dot(red_dichro,red_grating)),(np.dot(red_ccd,np.square(mirror)) * coating_eff)))
	red_noise = np.dot(counts_noise * red_total_eff_noise)

if (channel is 'blue') or (channel is 'both'):
	snr_blue = np.divide(blue_signal,np.sqrt(blue_signal + blue_noise + np.square(readnoise)))
if (channel is 'red') or (channel is 'both'):
	snr_red = np.divide(red_signal,np.sqrt(red_signal + red_noise + np.square(readnoise)))

if (channel is 'blue') or (channel is 'both'):
	sigma_blue = np.sqrt(blue_signal + blue_noise)
if (channel is 'red') or (channel is 'both'):
	sigma_red = np.sqrt(red_signal + red_noise)

if (channel is 'blue') or (channel is 'both'):
	error_blue = np.random.normal(loc=0, scale=(np.sqrt(blue_signal + blue_noise)),size=len(blue_snr))
if (channel is 'red') or (channel is 'both'):
	error_red = np.random.normal(loc=0, scale=(np.sqrt(red_signal + red_noise)),size=len(red_snr))


''' pre-plotting '''

thing_keys = ['snr','obv_spec','sky_background','dichroic_throughput','grating_throughput','ccd_qe','atmospheric_extinction']
plot_x = wavelength
if thing in thing_keys:
	if (thing.lower() is thing_keys[0]):
		title = 'Signal-to-Noise Ratio'
		labels = ['Angstrom','Signal-to-Noise Ratio']
		if (channel is 'blue') or (channel is 'both'):
			plot_y_blue = blue_snr
		if (channel is 'red') or (channel is 'both'):
			plot_y_red = red_snr
	elif (thing.lower() is thing_keys[1]):
		labels = ['Angstrom','Counts (#) per pixel']
		if noise:
			title = 'Observed Spectrum (with noise)'
			if (channel is 'blue') or (channel is 'both'):
				plot_y_blue = np.add(blue_signal,error_blue)
			if (channel is 'red') or (channel is 'both'):
				plot_y_red = np.add(red_signal,error_red)
		else:
			title = 'Observed Spectrum (without noise)'
			if (channel is 'blue') or (channel is 'both'):
				plot_y_blue = blue_signal
			if (channel is 'red') or (channel is 'both'):
				plot_y_red = red_signal
	elif (thing.lower() is thing_keys[2]):
		title = 'Observed Sky Background'
		labels = ['Angstrom','Counts (#) per pixel']
		if (channel is 'blue') or (channel is 'both'):
			plot_y_blue = blue_noise
		if (channel is 'red') or (channel is 'both'):
			plot_y_red = red_noise
	elif (thing.lower() is thing_keys[3]):
		title = 'Dichroic Throughput'
		labels = ['Angstrom','Throughput']
		if (channel is 'blue') or (channel is 'both'):
			plot_y_blue = blue_dichro
		if (channel is 'red') or (channel is 'both'):
			plot_y_red = red_dichro
	elif (thing.lower() is thing_keys[4]):
		title = 'Grating Throughput'
		labels = ['Angstrom','Throughput']
		if (channel is 'blue') or (channel is 'both'):
			plot_y_blue = blue_grating
		if (channel is 'red') or (channel is 'both'):
			plot_y_red = red_grating
	elif (thing.lower() is thing_keys[5]):
		title = 'CCD Quantum Efficiency'
		labels = ['Angstrom','QE']
		if (channel is 'blue') or (channel is 'both'):
			plot_y_blue = blue_ccd
		if (channel is 'red') or (channel is 'both'):
			plot_y_red = red_ccd
	elif (thing.lower() is thing_keys[6]):
		title = 'Atmospheric Extinction'
		labels = ['Angstrom','Throughput']
		plot_y_blue = spectres(wavelength,atmo_ext_x,atmo_ext_y)
		plot_y_red = plot_y_blue
	else:
		raise ValueError('{} Invalid thing ({})'.format(string_prefix,thing))
