#!/usr/bin/python3
'''
	GMACS ETC: ODC
		(omni-dimensional creator)
			builds the value volumes
				1) imports data
				2) reshapes data
				3) saves data

'''
import os
import time
import h5py
import numpy as np
import pandas as pd
from astropy import units as u
from tqdm import tqdm
from util import config as cfg

''' preamble '''
print('GMACS ETC: ODC')
time_start = time.time()

''' importing '''
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

# efficiency dots
grating_red_1 = np.dot(coalesced_grating_files[0]['a'],10)
grating_red_2 = coalesced_grating_files[0]['b']
grating_blue_1 = np.dot(coalesced_grating_files[1]['a'],10)
grating_blue_2 = coalesced_grating_files[1]['b']
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

''' dump re-shaped data to file '''
gdata = [coalesced_object_x,coalesced_object_y,coalesced_star_types,coalesced_filters,coalesced_sky_files,
		coalesced_grating_files,coalesced_ccd_files,coalesced_dichroic_files,coalesced_atmo_ext_files,grating_red_1,
		grating_red_2,grating_blue_1,grating_blue_2,ccd_efficiency_red_1,ccd_efficiency_red_2,ccd_efficiency_blue_1,
		ccd_efficiency_blue_2]	
fileout = h5py.File('data.h5','w')
fileout.create_dataset('etc_data',data=gdata)
fileout.close

''' afterword '''
print('\n[Success!] Build completed in {} seconds.'.format(str(time.time()-time_start)))