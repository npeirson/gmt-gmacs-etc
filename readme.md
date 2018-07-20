# GMACS-ETC 2.2
## Exposure Time Calculator for the Giant Magellan Telescope Multi-object Astronomical and Cosmological Spectrograph

Desiged for the Munnerlyn Astronomical Instrumentation Lab at Texas A&M University
Coded by the Kwisatz Haderach in 2018, based on previous versions by Ting Li

Documentation in `docs/`

TODOs:
- output all the re-parsed data so it doesn't have to be calculated every time
- add input disabling for plot calls that don't require them
- math problem somewhere in signal/noise ratio?
- add the loading screen (last)

Ok soldier listen up! Here's the query map!
```
		Tabs:
		 |-	SNR
		 |	 |- SNR_r
		 |	 |	 |- red_signal
		 |	 |	 |	 |- counts_red
		 |	 |	 |	 |	 |- power_red
		 |	 |	 |	 |	 |	 |- flux (function)
		 |	 |	 |	 |	 |	 |- area (user input)
		 |	 |	 |	 |	 |	 |- exp_time (user input)
		 |	 |	 |	 |	 |
		 |	 |	 |	 |	 |- wavelength (function)
		 |	 |	 |	 |	 |- constantants (h,c,1e10)
		 |	 |	 |	 |
		 |	 |	 |	 |- percent
		 |	 |	 |	 |	 |- funx
		 |	 |	 |	 |	 |	 |- sigma
		 |	 |	 |	 |	 |	 	 | - seeing (user input)
		 |	 |	 |	 |	 |
		 |	 |	 |	 |	 |- slit_size (user input)
		 |	 |	 |	 |	 |- seeing (user input)
		 |	 |	 |	 |
		 |	 |	 |	 |- red_total_eff
		 |	 |	 |		 |- red_dichro (file)
		 |	 |	 |		 |	 |- dichro_x
		 |	 |	 |		 |	 |- dichro_y2
		 |	 |	 |		 |
		 |	 |	 |		 |- grating_red
		 |	 |	 |		 |	 |- grating2 (file)
		 |	 |	 |		 |	 |- wavelength
		 |	 |	 |		 |- red_ccd
		 |	 |	 |		 |	 |- ccd2 (file)
		 |	 |	 |		 |	 |- wavelength
		 |	 |	 |		 |
		 |	 |	 |		 |- constants (coating eff)
		 |	 |	 |		 |- extinction (file)
		 |	 |	 |			 |- atmo_ext_x
		 |	 |	 |			 |- atmo_ext_y
		 |	 |	 |			 |- wavelength	 
		 |	 |	 |
		 |	 |	 |- red_noise
		 |	 |	 	 |- counts_noise_red
		 |	 |	 	 |	 |- sky_red
		 |	 |	 	 |	 	 |- sky_flux
		 |	 |	 	 |		 |	 |- sky_x (file)
		 |	 |	 	 |		 |	 |- sky_y (file)
		 |	 |	 	 |		 |	 |- wavelength
		 |	 |	 	 |		 |
		 |	 |	 	 |		 |- extension
		 |	 |	 	 |		 |	 |- seeing
		 |	 |	 	 |		 |	 |- slit_size
		 |	 |	 	 |		 |
		 |	 |	 	 |		 |- exp_time
		 |	 |	 	 |
		 |	 |	 	 |- red_total_eff_noise
		 |	 |			 |- red_dichro
		 |	 | 			 |- red_grating
		 |	 | 			 |- red_ccd
		 |	 | 			 |- coating_efficiency
		 |	 | 
		 |	 |- SNR_b
		 |	  	 |- same-same, but blue...
		 |	 	
		 |- Obs. Spectr. w/o Noise
		 |	 |- red_signal
		 |	 |- blue_signal
		 |	 |- wavelength
		 |	 
		 |-	Obs. Spectr. w/ Noise
		 |	 |- red_signal
		 |	 |- error_r
		 |	 |	 |- sigma_r
		 |	 |	 	 |- red_signal
		 |	 |	 	 |- red_noise
		 |	 |
		 |	 |- blue_signal
		 |	 |- error_b
		 |	 |- wavelength
		 |
		 |- Obs. Sky Background
		 |	 |- red_noise
		 |	 |- blue_noise
		 |	 |- wavelength
		 |
		 |-	Dichroic throughput
		 |	 |- red_dichro
		 |	 |	 |- dichro_x
		 |	 |	 |- dichro_y2
		 |	 |	 |- wavelength
		 |	 |
		 |	 |- blue_dichro
		 |	 |- wavelength
		 |
		 |-	Grating throughput
		 |	 |- red_grating
		 |	 |- blue_grating
		 |	 |- wavelength
		 |
		 |- CCD QE
		 |	 |- red_ccd
		 |	 |- blue_ccd
		 |	 |- wavelength
		 |
		 |-	Atmospheric Extinction
		 	 |- atmo_ext_x (file)
		 	 |- atmo_ext_y (file)
```
Clearly, then, wavelength is of the uppermost importance, followed by several other basic types.

TODO:
 [ ] update dichroic transition wavelength
 [ ] update telescope throughput
 [ ] add fiber transmission option
 [ ] update grating options and efficiencies
 [ ] update sky background wavelength coverage
 [ ] update optics throughput values
