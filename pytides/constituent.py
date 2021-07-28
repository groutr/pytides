
import string
import operator as op
import numpy as np
from . import nodal_corrections as nc
from functools import reduce

class BaseConstituent(object):
	xdo_int = {
		'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7, 'H': 8, 'I': 9,
		'J': 10, 'K': 11, 'L': 12, 'M': 13, 'N': 14, 'O': 15, 'P': 16, 'Q': 17,
		'R': -8, 'S': -7, 'T': -6, 'U': -5, 'V': -4, 'W': -3, 'X': -2, 'Y': -1,
		'Z': 0
	}

	int_xdo = {v:k for k, v in xdo_int.items()}

	def __init__(self, name, xdo='', coefficients=None, u=nc.u_zero, f=nc.f_unity):
		if coefficients is None:
			coefficients = []
		
		if xdo == '':
			self.coefficients = np.array(coefficients)
		else:
			self.coefficients = np.array(self.xdo_to_coefficients(xdo))
		self.name = name
		self.u = u
		self.f = f

	def xdo_to_coefficients(self, xdo):
		xd = [self.xdo_int[l.upper()] for l in xdo if l in string.ascii_letters]
		return xd

	def coefficients_to_xdo(self, coefficients):
		return ''.join(self.int_xdo[c] for c in coefficients)

	def V(self, astro):
		return np.dot(self.coefficients, self.astro_values(astro))

	def xdo(self):
		return self.coefficients_to_xdo(self.coefficients)

	def speed(self, a):
		return np.dot(self.coefficients, self.astro_speeds(a))

	def astro_xdo(self, a):
		return [a['T+h-s'], a['s'], a['h'], a['p'], a['N'], a['pp'], a['90']]

	def astro_speeds(self, a):
		return np.array([each.speed for each in self.astro_xdo(a)])

	def astro_values(self, a):
		return np.array([each.value for each in self.astro_xdo(a)])

	# Consider two out of phase constituents which travel at the same speed to
	# be identical
	def __eq__(self, c):
		return np.all(self.coefficients[:-1] == c.coefficients[:-1])

	def __hash__(self):
		return hash(tuple(self.coefficients[:-1]))

	def __str__(self):
		return self.name

	def __repr__(self):
		return f"BaseConstituent({self.name}, {self.u}, {self.f})"

class CompoundConstituent(BaseConstituent):

	def __init__(self, members=None, **kwargs):
		if members is None:
			self.members = []
		else:
			self.members = members

		if 'u' not in kwargs:
			kwargs['u'] = self.u
		if 'f' not in kwargs:
			kwargs['f'] = self.f

		super(CompoundConstituent,self).__init__(**kwargs)

		mcoefs = (c.coefficients * n for (c,n) in members)

		self.coefficients = reduce(op.add, mcoefs)

	def speed(self, a):
		return reduce(op.add, (n * c.speed(a) for (c,n) in self.members))

	def V(self, a):
		return reduce(op.add, (n * c.V(a) for (c,n) in self.members))

	def u(self, a):
		return reduce(op.add, (n * c.u(a) for (c,n) in self.members))

	def f(self, a):
		return reduce(op.mul, (c.f(a) ** abs(n) for (c,n) in self.members))

	def __repr__(self):
		return f"CompoundConstituent({self.name}, {self.members})"


# # # # # #  Base Constituents
# Long Term
_Z0 = BaseConstituent(name = 'Z0', xdo = 'Z ZZZ ZZZ' )
_Sa = BaseConstituent(name = 'Sa', xdo = 'Z ZAZ ZZZ' )
_Ssa = BaseConstituent(name = 'Ssa', xdo = 'Z ZBZ ZZZ')
_Mm = BaseConstituent(name = 'Mm', xdo = 'Z AZY ZZZ', f = nc.f_Mm)
_Mf = BaseConstituent(name = 'Mf', xdo = 'Z BZZ ZZZ', u = nc.u_Mf, f = nc.f_Mf)

# Diurnals
_Q1 = BaseConstituent(name = 'Q1', xdo = 'A XZA ZZA', u = nc.u_O1, f = nc.f_O1)
_O1 = BaseConstituent(name = 'O1', xdo = 'A YZZ ZZA', u = nc.u_O1, f = nc.f_O1)
_K1 = BaseConstituent(name = 'K1', xdo = 'A AZZ ZZY', u = nc.u_K1, f = nc.f_K1)
_J1 = BaseConstituent(name = 'J1', xdo = 'A BZY ZZY', u = nc.u_J1, f = nc.f_J1) #?

# M1 is a tricky business for reasons of convention, rather than theory.  The
# reasons for this are best summarised by Schureman paragraphs 126, 127 and in
# the comments found in congen_input.txt of xtides, so I won't go over all this
# again here.

_M1 = BaseConstituent(name = 'M1', xdo = 'A ZZZ ZZA', u = nc.u_M1, f = nc.f_M1)
_P1 = BaseConstituent(name = 'P1', xdo = 'A AXZ ZZA')
_S1 = BaseConstituent(name = 'S1', xdo = 'A AYZ ZZZ')
_OO1 = BaseConstituent(name = 'OO1', xdo = 'A CZZ ZZY', u = nc.u_OO1, f = nc.f_OO1)

# Semi-Diurnals
_2N2 = BaseConstituent(name = '2N2', xdo = 'B XZB ZZZ', u = nc.u_M2, f = nc.f_M2)
_N2 = BaseConstituent(name = 'N2', xdo = 'B YZA ZZZ', u = nc.u_M2, f = nc.f_M2)
_nu2 = BaseConstituent(name = 'nu2', xdo = 'B YBY ZZZ', u = nc.u_M2, f = nc.f_M2)
_M2 = BaseConstituent(name = 'M2', xdo = 'B ZZZ ZZZ', u = nc.u_M2, f = nc.f_M2)
_lambda2 = BaseConstituent(name = 'lambda2', xdo = 'B AXA ZZB', u = nc.u_M2, f = nc.f_M2)
_L2 = BaseConstituent(name = 'L2', xdo = 'B AZY ZZB', u = nc.u_L2, f = nc.f_L2)
_T2 = BaseConstituent(name = 'T2', xdo = 'B BWZ ZAZ')
_S2 = BaseConstituent(name = 'S2', xdo = 'B BXZ ZZZ')
_R2 = BaseConstituent(name = 'R2', xdo = 'B BYZ ZYB')
_K2 = BaseConstituent(name = 'K2', xdo = 'B BZZ ZZZ', u = nc.u_K2, f = nc.f_K2)

# Third-Diurnals
ufunc = lambda a:nc.u_Modd(a,3)
ffunc = lambda a:nc.f_Modd(a,3)
_M3 = BaseConstituent(name = 'M3', xdo = 'C ZZZ ZZZ', u = ufunc, f = ffunc)


# # # # # #  Compound Constituents
# Long Term
_MSF = CompoundConstituent(name = 'MSF', members = [(_S2, 1), (_M2, -1)])
_MSTM = CompoundConstituent(name = 'MSTM', members = [(_Mf, 2), (_Mm, -1)])

# Diurnal
_2Q1 = CompoundConstituent(name = '2Q1', members = [(_N2, 1), (_J1, -1)])
_rho1 = CompoundConstituent(name = 'rho1', members = [(_nu2, 1), (_K1, -1)])

# Semi-Diurnal
_mu2 = CompoundConstituent(name = 'mu2', members = [(_M2, 2), (_S2, -1)]) # 2MS2
_2SM2 = CompoundConstituent(name = '2SM2', members = [(_S2, 2), (_M2, -1)])
_OP2 = CompoundConstituent(name = 'OP2', members = [(_P1, 1), (_O1, 	1)])
_MNS2 = CompoundConstituent(name = 'MNS2', members = [(_S2, -1), (_M2, 1), (_N2, 1)])

# Third-Diurnal
_2MK3 = CompoundConstituent(name = '2MK3', members = [(_M2, 1), (_O1, 1)])
_MK3 = CompoundConstituent(name = 'MK3', members = [(_M2, 1), (_K1, 1)])
_MO3 = CompoundConstituent(name = 'MO3',  members = [(_M2, 1), (_O1, 1)]) #(?)

# Quarter-Diurnal
_MN4 = CompoundConstituent(name = 'MN4', members = [(_M2, 1), (_N2, 1)])
_M4 = CompoundConstituent(name = 'M4', members = [(_M2, 2)])
_MS4 = CompoundConstituent(name = 'MS4', members = [(_M2, 1), (_S2, 1)])
_S4 = CompoundConstituent(name = 'S4', members = [(_S2, 2)])
_N4 = CompoundConstituent(name = 'N4', members = [(_N2, 2)])
_L4 = CompoundConstituent(name = 'L4', members = [(_L2, 2)])
_SN4 = CompoundConstituent(name = 'SN4', members = [(_S2, 1), (_N2, 1)])
_MK4 = CompoundConstituent(name = 'MK4', members = [(_M2, 1), (_K2, 1)])
_2MSN4 = CompoundConstituent(name = '2MSN4', members = [(_S2, 1), (_M2, 2), (_N2, -1)])

# Sixth-Diurnal
_M6 = CompoundConstituent(name = 'M6', members = [(_M2, 3)])
_S6 = CompoundConstituent(name = 'S6', members = [(_S2, 3)])

# Eighth-Diurnals
_M8 = CompoundConstituent(name = 'M8', members = [(_M2, 4)])
_S8 = CompoundConstituent(name = 'S8', members = [(_S2, 4)])

# Tenth-Diurnals
_M10 = CompoundConstituent(name = 'M10', members = [(_M2, 5)])

# Twelveth-Diurnals
_M12 = CompoundConstituent(name = 'M12', members = [(_M2, 6)])


noaa = [
	_M2, _S2, _N2, _K1, _M4, _O1, _M6, _MK3, _S4, _MN4, _nu2, _S6, _mu2, _2N2,
	_OO1, _lambda2, _S1, _M1, _J1, _Mm, _Ssa, _Sa, _MSF, _Mf,	_rho1, _Q1, _T2,
	_R2, _2Q1, _P1, _2SM2, _M3, _L2, _2MK3, _K2 ,_M8, _MS4, _MO3, _SN4, _MNS2,
	_N4, _L4, _M10, _M12, _S8, _MK4, _OP2,_2MSN4, _MSTM,
]

_constit_mapping = {
	"Z0": _Z0,
	"SA": _Sa,
	"SSA": _Ssa,
	'MM': _Mm,
	"MF": _Mf,
	"Q1": _Q1,
	"O1": _O1,
	"K1": _K1,
	"J1": _J1,
	"M1": _M1,
	"P1": _P1,
	"S1": _S1,
	"OO1": _OO1,
	"2N2": _2N2,
	"N2": _N2,
	"NU2": _nu2,
	"M2": _M2,
	"LAMBDA2": _lambda2,
	"L2": _L2,
	"T2": _T2,
	"S2": _S2,
	"R2": _R2,
	"K2": _K2,
	"M3": _M3,
	"MSF": _MSF,
	"MSTM": _MSTM,
	"2Q1": _2Q1,
	"RHO1": _rho1,
	"MU2": _mu2,
	"2SM2": _2SM2,
	"OP2": _OP2,
	"MNS2": _MNS2,
	"2MK3": _2MK3,
	"MK3": _MK3,
	"MO3": _MO3,
	"MN4": _MN4,
	"M4": _M4,
	"MS4": _MS4,
	"S4": _S4,
	"N4": _N4,
	"L4": _L4,
	"SN4": _SN4,
	"MK4": _MK4,
	"2MSN4": _2MSN4,
	"M6": _M6,
	"S6": _S6,
	"M8": _M8,
	"S8": _S8,
	"M10": _M10,
	"M12": _M12
}

def get_constituent_by_name(constituent_name):
	## Upper-case the provided name to match the mapping keys
	cn = constituent_name.upper()
	if cn in _constit_mapping.keys():
		return _constit_mapping[cn]
	
	## Not found, raise error
	raise ValueError("No constituent found for name '{}'".format(constituent_name))
	

def get_constituent_names():
	return list(_constit_mapping.keys())
	

#
#  u'RO1', u'SIGMA1', u'PI1',
# u'FI1', u'CHI1', u'THETA1', u'SO1', u'MP1', u'PSI1',
# , u'KJ2', u'MSN2', u'MKS2', u'OQ2',
# u'SK3', u'SO3', u'SK4', u'2MS6', u'2MN6',
# u'2SM6', u'MSN6', u'2MK6', u'MSK6', u'SK2', u'S3',
#  u'MSM', u'MTM', u'MSQM', u'M(SK)2', u'M(KS)2',
# u'3MS8', u'2(MS)8', u'2MSN8', u'2MNS4', u'3MK4', u'3MS4', u'2MSK4',
# u'2MKS4', u'3MN4', u'2SMK4', u'2MSN4', u'MNK6', u'MSNK8', u'2SNM4',
# u'SP3', u'KQ1', u'MQ3', u'K3', u'3SM4', u'2SN2', u'2SN6',
# , u'NO1', u'2MS2', u'ALP1', u'EPS2', u'ETA2', u'UPS1',
# u'2MN2', u'3MK7', u'2MK5', u'2SK5', u'TAU1', u'H1', u'H2', u'4MS6',
# u'MQM', u'SIG1', u'M11', u'KI1', u'PHI1', u'TTA1', u'E2',
# u'2MK2', u'MSK2', u'LA2', u'A0', u'3MN8', u'2MSK8'
