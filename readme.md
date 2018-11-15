# GMACS ETC

server version just needs toggle channel option to get the glyphs to initialize, and there are some bugs in various buttons and such. Also, the layout is pretty messy too, so it isn't just showing up poorly in your browser; it needs some changes. I think they're because I haven't moved over the script that manages what the dropdown menus display and such. That's what the (currently empty) `enabler` function is for.

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

