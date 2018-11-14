#!/usr/bin/python3
"""
	etc tools
	----------------
	re-usable and unwieldly functions that are useful,
	and also best kept in isolation.
	Like mushrooms, poets, and fuel rods.

"""

import numpy as np

import keys as etkeys
import defaults as dfs
'''
def registrar(kwargs):
	validated_args,validated_values = np.empty(0),np.empty(0)
	for _key,_value in kwargs.items():
		if isinstance(_key,str):
			_key,_value = _key.lower(),_value.lower()
			for i in range(len(etkeys.arguments)): 		# go through list of lists
				if _key in etkeys.arguments[i]:			# see if key is in each one
					_key = etkeys.arguments[i][0]		# if it is, standardize it
				else:
					raise ValueError("{} Registrar reports key not in keychain!".format(dfs.string_prefix))

		if isinstance(_value,str):
			if (_key == 0): # mode
				for j in range(len(etkeys.mode)):
					if _value in etkeys.mode[j]:
						_value = int(np.where(np.asarray(etkeys.mode[j])==_value)[0])
					else:
						raise ValueError("{} Registrar reports value not in keychain!".format(dfs.string_prefix))
'''







def value_registrar(_obj):
	if isinstance(_obj,str):
		_obj = _obj.lower()
	if (_obj in np.asarray(etkeys.keychain).flatten()) or ((isinstance(_obj,float)) and ((_obj >= 0) and(_obj <= 13000.))):
		for i in range(len(etkeys.keychain)):
			if _obj in np.asarray(etkeys.keychain[i]).flatten():
				for j in range(len(etkeys.keychain[i])):
					if (_obj == _inst) or (_obj == j) or (isinstance(_obj,float)):
						return j
	else:
		raise ValueError("{} Registrar reports value not in keychain!".format(dfs.string_prefix))


def arg_registrar(_arg):
	if isinstance(_arg,str):
		_arg = _arg.lower()
	if _arg in np.asarray(etkeys.arguments).flatten():
		for i,argument_list in enumerate(etkeys.arguments):
			if _arg in argument_list:
				return etkeys.arguments[i][0]
	else:
		raise ValueError("{} Registrar reports argument not in keychain!".format(dfs.string_prefix))


def validate_args(kwargs):
	validated_args,validated_values = np.empty(0),np.empty(0)
	for _key,_value in kwargs.items():
		if isinstance(_value,bool):
			if _key in etkeys.arguments[-1]:
				validated_args = np.append(validated_args,_key)
				validated_values = np.append(validated_values,_value)
			else:
				raise ValueError("{} Subsiste sermonem statim (sss) must be type `bool`.".format(dfs.string_prefix))
		elif isinstance(_value,list) or isinstance(_value,np.ndarray) or isinstance(_value,range):
			if _key in etkeys.arguments[2]: # wavelength
				if (_value[0] >= dfs.wavelength_limits[0]) and (_value[-1] <= dfs.wavelength_limits[1]):
					valid_wavelength = np.array(_value)
					#validated_args = np.append(validated_args,_key)
					#validated_values = np.append(validated_values,[np.array(_value)])
				else:
					raise ValueError("{} Requested wavelength range exceeds calculable range ({}--{} nm)".format(dfs.string_prefix,dfs.wavelength_limits[0],dfs.wavelength_limits[1]))
			else:
				raise TypeError("{} Only attribute `wavelength` may be array-like.".format(dfs.string_prefix))
		elif isinstance(_value,int) or isinstance(_value,float):
			if _key in etkeys.arguments[3]: # exposure time
				if (_value > 0.) and (_value <= 2.628e6): # limit to one month...
					validated_args = np.append(validated_args,_key)
					validated_values = np.append(validated_values,_value)
				else:
					raise ValueError("{} Requested exposure time ({} sec) exceeds limit (2.628e6 sec)".format(dfs.string_prefix,_value))
			elif _key in etkeys.arguments[7]: # magntidue
				if (_value >= -50.) and (_value <= 50.):
					validated_args = np.append(validated_args,_key)
					validated_values = np.append(validated_values,_value)
				else:
					raise ValueError("{} Requested magnitude ({}) is outside acceptable range (+/- 50)".format(dfs.string_prefix,_value))
			elif _key in etkeys.arguments[8]: # redshift
				if (_value >= 0.) and (_value <= 1090.):
					validated_args = np.append(validated_args,_key)
					validated_values = np.append(validated_values,_value)
				else:
					raise ValueError("{} Requested redshift ({}) is outside acceptable range (0--12)".format(dfs.string_prefix,_value))
			elif _key in etkeys.arguments[9]: # seeing
				if (_value >= 0) and (_value <= 50):
					validated_args = np.append(validated_args,_key)
					validated_values = np.append(validated_values,_value)
				else:
					raise ValueError("{} Requested seeing ({}) is outside acceptable range (0--50 arcsec)".format(dfs.string_prefix,_value))
			elif _key in etkeys.arguments[10]: # slit width
				validated_args = np.append(validated_args,_key)
				validated_values = np.append(validated_values,_value)
			elif _key in etkeys.arguments[11]: # moon days
				if int(_value) in etkeys.moon_days[0]:
					validated_args = np.append(validated_args,_key)
					validated_values = np.append(validated_values,_value)
				else:
					raise ValueError("{} Requested `days since new moon` is not a valid number. Acceptable options are: {}".format(dfs.string_prefix,etkeys.moon_days))
			elif _key in etkeys.arguments[14]: # bin_opt
				if int(_value) in etkeys.bin_opt[0]:
					validated_args = np.append(validated_args,_key)
					validated_values[int(np.where(np.asarray(validated_args)==_key)[0])] = _value # move this
				else:
					raise ValueError("{} Requested pixel binning invalid! Please select from: {}".format(dfs.string_prefix,etkeys.bin_opt))
			else:
				pass # ok nevermind then

		elif isinstance(_value,str) or (_value == None): # catch all the rest
			if _key in etkeys.arguments[4]: # object_type
				for i in range(len(etkeys.object_type)):
					_value = (i,np.where(np.asarray(etkeys.object_type[i])==_value))

			for i in range(len(etkeys.arguments)):
				if _key in etkeys.arguments[i]:
					_key = etkeys.arguments[i][0]
					validated_args = np.append(validated_args,_key)

			for j in range(len(etkeys.keychain)):		
				for k in range(len(etkeys.keychain[j])):
					if _value in etkeys.keychain[j][k]:
						if (len(etkeys.keychain[j][k]) > 1):
							for l in range(len(etkeys.keychain[j][k])):
								if (etkeys.keychain[j][k][l] == _value):
									print('replacing {} ith {}'.format(_value,l))
									validated_values = np.append(validated_values,l)
									break
						elif (len(etkeys.keychain[j][k]) == 1):
							if (_value == etkeys.keychain[j][0][0]):
								break
						else:
							raise ValueError("{} Invalid value: ({})".format(dfs.string_prefix,_value))
		else:
			raise ValueError("{} Invalid argument: ({})".format(dfs.string_prefix,_key))	
	for i in range(len(validated_args)):
		print(validated_args[i])
	for j in range(len(validated_values)):
		print(validated_values[j])
	#print("arg: {}\tval: {}".format(validated_args[i],validated_values[i])) # for debugging

	# final lenegth check	
	if (len(validated_args) != len(validated_values)):
		#print("{} : {}".format(validated_args,validated_values))
		raise ValueError("{} Argument / value mis-match! {} != {}".format(dfs.string_prefix,len(validated_args),len(validated_values)))
	else:
		#for i in range(len(validated_args)): # print("arg: {}\tval: {}".format(validated_args[i],validated_values[i])) # for debugging			
		return dict(zip(validated_args,validated_values)), valid_wavelength


def validate_points(x,y):
	if (isinstance(x,np.ndarray) and isinstance(y,np.ndarray)) or (isinstance(np.asarray(x),np.ndarray) and isinstance(np.asarray(y),np.ndarray)):
		if (x.shape[0] == y.shape):
			return x,y
		elif (x.shape[0] == (y.shape[0] / 2)):
			return x,y[0],y[1]
		else:
			raise ValueError("{} Plot columns must be equal in length. (x={},y={})".format(dfs.string_prefix,x.shape[0],y.shape[0]))
	else:
		raise TypeError("{} Invalid points proposed for plotting! Must be type `nd.array`".format(dfs.string_prefix))
