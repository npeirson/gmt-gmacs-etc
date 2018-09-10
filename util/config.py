#!/usr/bin/python3

'''
	GMACS ETC : Congifuration file

'''


''' paths '''
# explicit because have to be sure indices are consistent

galaxy_types_path = 'core/kinney/'
galaxy_types_files = ['SB1.txt','SB2.txt','SB3.txt','SB4.txt','SB5.txt','SB6.txt','S0.txt','Sa.txt',
						'Sb.txt','Sc.txt','bulge.txt','ellipticals.txt','lbg_all_flam.txt']
star_types_path = 'core/pickle/'
star_types_files = ['o5v.dat','b0v.dat','b57v.dat','a0v.dat','a5v.dat',
					'f0v.dat','g0v.dat','g5v.dat','k0v.dat','k5v.dat',
					'm0v.dat','m5v.dat']
filter_path = 'core/filters/'
filter_files = ['photonUX.txt','photonB.txt','photonV.txt','photonR.txt','photonI.txt','u.txt','g.txt','r.txt','i.txt','z.txt']

skyfiles_path = 'core/skybackground/'
skyfiles_files = ['00d.txt','03d.txt','07d.txt','10d.txt','14d.txt']

efficiency_path = 'core/efficiencies/'
efficiency_grating_files = ['grating_red_low_res.txt','grating_blue_low_res.txt']
efficiency_ccd_files = ['e2v_red.txt','e2v_blue.txt']
dichroic_path = 'core/'
dichroic_files = ['dichroic.txt']
atmo_ext_path = 'core/'
atmo_ext_files = ['atmo_extinction.dat']