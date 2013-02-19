import scipy
import constants

## constants
#coulombPerElectron = 1.60217733e-19 
#nm = 1e-9

# simulation setup
coulombPotetial = 1

## careful when using 2D as it is not debugged yet
#dim = 3

#particles = 0

x = 2500.*constants.nm
y = x
z = 250.*constants.nm
r = x #parameter used to determine the radial confinement potential

periodic = scipy.array([0., 0., 0.])

#shapeDeposit = 'rectangel'
#xDeposit = x/2.
#yDeposit = y/2.
zDeposit = z

shapeDeposit = 'circle'
rDeposit = x/10.

#shapeDeposit = 'repeat'

#V_z = .1
V_x_constant = 0.
V_y_constant = 0.
V_z_constant = .25
# roughly 1kV per 2.54cm
#V_z = 0.02
V_x_linear = 0.
V_y_linear = 0.
V_z_linear = 0.
V_r_linear = 10.

# Field over the full box
E_x_constant = V_x_constant/(2.*x)
E_y_constant = V_y_constant/(2.*y)
E_z_constant = V_z_constant/(2.*z)
# Field over half the box
E_x_linear = 2.*V_x_linear/(x**2)
E_y_linear = 2.*V_y_linear/(y**2)
E_z_linear = 2.*V_z_linear/(z**2)
E_r_linear = 2.*V_r_linear/(r**2)

CosntantE = scipy.array([E_x_constant, E_y_constant, E_z_constant])
LinearE = scipy.array([E_x_linear, E_y_linear, E_z_linear])
RadialLinearE = E_r_linear

iterations = 100
sampling = 10

timeStep = 2.5e-16
timeDeposite = 400.*timeStep
#ensure that a particle initially is deposited
timeLastDeposit = 2.*timeDeposite
