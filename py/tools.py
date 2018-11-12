
import numpy as np

import keys as etkeys
import defaults as dfs


def value_registrar(_obj):
	if isinstance(_obj,str):
		_obj = obj.lower()

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
		_arg = arg.lower()

	if _arg in np.asarray(etkeys.arguments).flatten():
		for i,argument_list in enumerate(etkeys.arguments):
			if _arg in argument_list:
				return etkeys.arguments[i][0]
	else:
		raise ValueError("{} Registrar reports argument not in keychain!".format(dfs.string_prefix))


