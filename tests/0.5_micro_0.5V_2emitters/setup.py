import scipy
import constants

sampling = 1000
timeStep = 2.5e-16

# simulation setup
coulombInteractions = 1

###############
## BOX START ##
x = 2500.*constants.nm
y = x
z = 250.*constants.nm # half box length
r = x # parameter used to determine the radial confinement potential
periodic = scipy.array([0., 0., 0.])
##  BOX END  ##
###############

####################
## EMISSION START ##
emissionNrAttempts = 100
emissionType = 'space_charge'
steps = 20000
emissionExtraSteps = 0
emissionShape = 'circle'
emissionR = x/50.
emissionZ = z
##  EMISSION END  ##
####################

##################################
## CONFINEMENT POTENTIALS START ##
# constant terms
V_xConstant = 0.
V_yConstant = 0.
V_zConstant = .5
# linear terms
V_xLinear = 0.
V_yLinear = 0.
V_zLinear = 0.
#V_rLinear = 10.
V_rLinear = 0.

# corresponding fields computed
# constant terms
E_xConstant = V_xConstant/(2.*x)
E_yConstant = V_yConstant/(2.*y)
E_zConstant = V_zConstant/(2.*z)
# linear terms
E_xLinear = 2.*V_xLinear/(x**2)
E_yLinear = 2.*V_yLinear/(y**2)
E_zLinear = 2.*V_zLinear/(z**2)
E_rLinear = 2.*V_rLinear/(r**2)
# fields collected
constantE = scipy.array([E_xConstant, E_yConstant, E_zConstant])
linearE = scipy.array([E_xLinear, E_yLinear, E_zLinear])
radialLinearE = E_rLinear
##  CONFINEMENT POTENTIALS END  ##
##################################

