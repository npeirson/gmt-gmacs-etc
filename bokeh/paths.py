#!/usr/bin/python3

import os

''' paths '''
# get path to values.py
absFilePath = os.path.abspath(__file__)
# get directory
fileDir = os.path.dirname(os.path.abspath(__file__))
# get parent directory
parentDir = os.path.dirname(fileDir)

galaxy_path = parentDir + '/core/kinney/'
stellar_path = parentDir + '/core/pickle/'
skyfiles_path = parentDir + '/core/skybackground/'
filter_path = parentDir + '/core/filters/'
atmo_ext_path = parentDir + '/core/atmo_extinction.dat'
dichroic_path = parentDir + '/core/dichroic5577.txt'
grating_path = [parentDir + '/core/grating_red_low_315-1100.txt',parentDir + '/core/grating_blue_low_315-1100.txt']
ccd_path = [parentDir + '/core/e2v-astro-multi-2-DD.txt',parentDir + '/core/e2v_blue.txt']
vega_file = parentDir + '/core/alpha_lyr_stis_005.ascii'
mirror_file = parentDir + '/core/mirror_300-1200nm.txt'
galaxy_files = ['SB1','SB2','SB3','SB4','SB5','SB6','S0','Sa','Sb','Sc','bulge','ellipticals','lbg_all_flam']
stellar_files = ['o5v.dat','b0v.dat','b57v.dat','a0v.dat','a5v.dat','f0v.dat','g0v.dat','g5v.dat','k0v.dat','k5v.dat','m0v.dat','m5v.dat','f5v.dat']
skyfiles = ['00d_315-1200nm.csv','03d_315-1200nm.csv','07d_315-1200nm.csv','10d_315-1200nm.csv','14d_315-1200nm.csv']
filter_files = ['BesB.dat','BesI.dat','BesR.dat','BesU.dat','BesV.dat','g.dat','i.dat','r.dat','u.dat','z.dat','photonb.dat','photonI.dat','photonR.dat','photonUX.dat','photonV.dat']
#filter_files = ['photonUX.dat','photonB.dat','photonV.dat','photonR.dat','photonI.dat','u.dat','g.dat','r.dat','i.dat','z.dat']