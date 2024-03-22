
import math

# The following functions take a dictionary of astronomical values (degrees)
# and return dimensionless scale factors for constituent amplitudes.

def f_unity(a):
	return 1.0

# Schureman equations 73, 65
def f_Mm(a):
	omega = math.radians(a['omega'].value)
	i = math.radians(a['i'].value)
	I = math.radians(a['I'].value)
	mean = (2/3.0 - math.sin(omega) ** 2)*(1 - 3/2.0 * math.sin(i) ** 2)
	return (2/3.0 - math.sin(I)**2) / mean

# Schureman equations 74, 66
def f_Mf(a):
	omega = math.radians(a['omega'].value)
	i = math.radians(a['i'].value)
	I = math.radians(a['I'].value)
	mean = math.sin(omega)**2 * math.cos(0.5*i)**4
	return math.sin(I)**2 / mean

# Schureman equations 75, 67
def f_O1(a):
	omega = math.radians(a['omega'].value)
	i = math.radians(a['i'].value)
	I = math.radians(a['I'].value)
	mean = math.sin(omega) * math.cos(0.5*omega)**2 * math.cos(0.5*i)**4
	return (math.sin(I) * math.cos(0.5*I)**2) / mean

# Schureman equations 76, 68
def f_J1(a):
	omega = math.radians(a['omega'].value)
	i = math.radians(a['i'].value)
	I = math.radians(a['I'].value)
	mean = math.sin(2*omega) * (1-3/2.0 * math.sin(i)**2)
	return math.sin(2*I) / mean

# Schureman equations 77, 69
def f_OO1(a):
	omega = math.radians(a['omega'].value)
	i = math.radians(a['i'].value)
	I = math.radians(a['I'].value)
	mean = math.sin(omega) * math.sin(0.5*omega)**2 * math.cos(0.5*i)**4
	return math.sin(I) * math.sin(0.5*I)**2 / mean

# Schureman equations 78, 70
def f_M2(a):
	omega = math.radians(a['omega'].value)
	i = math.radians(a['i'].value)
	I = math.radians(a['I'].value)
	mean = math.cos(0.5*omega)**4 * math.cos(0.5*i)**4
	return math.cos(0.5*I)**4 / mean

# Schureman equations 227, 226, 68
# Should probably eventually include
# the derivations of the magic numbers (0.5023 etc).
def f_K1(a):
	omega = math.radians(a['omega'].value)
	i = math.radians(a['i'].value)
	I = math.radians(a['I'].value)
	nu = math.radians(a['nu'].value)
	sin2Icosnu_mean = math.sin(2*omega) * (1-3/2.0 * math.sin(i)**2)
	mean = 0.5023*sin2Icosnu_mean + 0.1681
	return math.sqrt(0.2523*math.sin(2*I)**2 +
			0.1689*math.sin(2*I)*math.cos(nu)+0.0283) / mean

# Schureman equations 215, 213, 204
# It can be (and has been) confirmed that the exponent for R_a reads 1/2
#  via Schureman Table 7
def f_L2(a):
	P = math.radians(a['P'].value)
	I = math.radians(a['I'].value)
	_tanI2 = math.tan(0.5 * I)
	R_a_inv = math.sqrt(1 -
			12*_tanI2**2 * math.cos(2*P)+
			36*_tanI2**4)
	return f_M2(a) * R_a_inv

# Schureman equations 235, 234, 71
# Again, magic numbers
def f_K2(a):
	omega = math.radians(a['omega'].value)
	i = math.radians(a['i'].value)
	I = math.radians(a['I'].value)
	_sinI = math.sin(I)
	nu = math.radians(a['nu'].value)
	sinsqIcos2nu_mean = math.sin(omega)**2 * (1-3/2.0 * math.sin(i)**2)
	mean = 0.5023*sinsqIcos2nu_mean + 0.0365
	return math.sqrt(0.2523*_sinI**4 +
	 		0.0367* _sinI**2 *math.cos(2*nu)+0.0013) / mean

# Schureman equations 206, 207, 195
def f_M1(a):
	P = math.radians(a['P'].value)
	I = math.radians(a['I'].value)
	_cosI = math.cos(I)
	_cosI2 = math.cos(0.5 * I)
	Q_a_inv = math.sqrt(0.25 +
				1.5*_cosI*math.cos(2*P)*_cosI2**(-0.5) +
				2.25*_cosI**2 * _cosI2**(-4))
	return f_O1(a) * Q_a_inv

# See e.g. Schureman equation 149
def f_Modd(a, n):
	return f_M2(a) ** (n / 2.0)

# Node factors u, see Table 2 of Schureman.

def u_zero(a):
	return 0.0

def u_Mf(a):
	return -2.0 * a['xi'].value

def u_O1(a):
	return 2.0 * a['xi'].value - a['nu'].value

def u_J1(a):
	return -a['nu'].value

def u_OO1(a):
	return -2.0 * a['xi'].value - a['nu'].value

def u_M2(a):
	return 2.0 * a['xi'].value - 2.0 * a['nu'].value

def u_K1(a):
	return -a['nup'].value

# Schureman 214
def u_L2(a):
	I = math.radians(a['I'].value)
	P = math.radians(a['P'].value)
	_2P = 2 * P
	R = math.degrees(math.atan(
		math.sin(_2P)/(1/6.0 * math.tan(0.5*I) **(-2) -math.cos(_2P))
		))
	return 2.0 * a['xi'].value - 2.0 * a['nu'].value - R

def u_K2(a):
	return -2.0 * a['nupp'].value

# Schureman 202
def u_M1(a):
	I = math.radians(a['I'].value)
	P = math.radians(a['P'].value)
	_cosI = math.cos(I)
	Q = math.degrees(math.atan((5*_cosI-1)/(7*_cosI+1)*math.tan(P)))
	return a['xi'].value - a['nu'].value + Q

def u_Modd(a, n):
	return n/2.0 * u_M2(a)
