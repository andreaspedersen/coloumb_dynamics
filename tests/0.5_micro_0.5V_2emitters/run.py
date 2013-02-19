import scipy
from scipy import stats
import pylab
import random
import pickle

import particles
reload(particles)

import movers
reload(movers)

import emitters
reload(emitters)

import boundaries
reload(boundaries)

import helper_functions
reload(helper_functions)

import constants
reload(constants)

import plot_tools
reload(plot_tools)

import analyze_tools
reload(analyze_tools)

import setup
reload(setup)


#----=======#######SETUP#######=======----
p = particles.Particles(0, setup.coulombInteractions)
p.set_box(scipy.array([setup.x, setup.y, setup.z])*setup.periodic)

activeEmitters = []
activeBoundaries = []

#----=======#######SETUP EXTERNAL POTENTIAL#######=======----
E_z = setup.V_zConstant/(2.*setup.z)
CosntantE = scipy.array([0., 0., E_z])

if setup.constantE.any():
    p.set_external_potential(setup.constantE, "Constant")

if setup.linearE.any():
    p.set_external_potential(setup.linearE, "Linear")

if setup.radialLinearE:
    p.set_external_potential(scipy.array([setup.radialLinearE, 0., 0.]), "RadialLinear")

#----=======#######EMITTER#######=======----	
eZ = setup.emissionZ
	
eR = setup.emissionR

active = {}
active['emitter'] = emitters.PlaneSpaceCharge(-eZ, eR, -500.*constants.nm, 0.*constants.nm)
active['delay'] = 0
activeEmitters.append(active)
active = {}
active['emitter'] = emitters.PlaneSpaceCharge(-eZ, eR, 500.*constants.nm, 0.*constants.nm)
active['delay'] = 10
activeEmitters.append(active)
		
#----=======#######SETUP BOUNDARIES#######=======----
x = setup.x
y = setup.y
z = setup.z

activeBoundaries.append(boundaries.Min('Reflecting', 'x', -x))
activeBoundaries.append(boundaries.Max('Reflecting', 'x',  x))
activeBoundaries.append(boundaries.Min('Reflecting', 'y', -y))
activeBoundaries.append(boundaries.Max('Reflecting', 'y',  y))
activeBoundaries.append(boundaries.Min('Reflecting', 'z', -z))
activeBoundaries.append(boundaries.Max('Absorbing',  'z',  z))

#----=======#######SETUP DRIVING CODE#######=======----
dataAll = []
posAll = []

phase = None
if phase is not None:
    phase['r'] = []
    phase['v'] = []
    phase['a'] = []

timeLastDeposit = 0
steps = 0

allAbsorbed = {}
allAbsorbed['nr'] = []
allAbsorbed['coordinates'] = []

allEmitted = {}
allEmitted['nr'] = []
allEmitted['coordinates'] = []

nr = []
totalTime = list([0.0])

m = movers.VelocityVerlet(p, setup.timeStep)
data = [timeLastDeposit, steps, posAll, allAbsorbed, allEmitted]


#----=======#######EMIT ELECTRONS#######=======----	
def emit(data):
    timeLastDeposit = data[0]
    steps = data[1]
    posAll = data[2]
    allAbsorbed = data[3]
    allEmitted = data[4]

    extraSteps = setup.emissionExtraSteps
    moreIterations = 1

    while moreIterations:
        # impose boundary conditions
        for boundary in activeBoundaries:
            absorbedStep = boundary.check_all(p)
            # keep track on how many atoms that have been absorbed
            if absorbedStep:
                for absorbed in absorbedStep:
                    allAbsorbed['coordinates'].append(absorbed)
        allAbsorbed['nr'].append(len(allAbsorbed['coordinates']))

	# emit electrons
        allDone = 0
        for emitter in activeEmitters:
            emitter['active'] = 1
 
        while not allDone:
	    allDone = 1
            for active in activeEmitters:
                if active['active'] and active['delay'] < steps:
                    emitted = active['emitter'].add_particle(p)
                    if emitted == None:
                        active['active'] = 0
                    else:
                        allDone = 0
                        allEmitted['coordinates'].append(emitted[0])


        allEmitted['nr'].append(len(allEmitted['coordinates']))
 
        totalTime.append(totalTime[-1] + setup.timeStep)
                    
        if steps%setup.sampling == 0:
            nr.append(len(p))
            if posAll is not None:
                pos = p.get_R()
                posAll.append(pos)
            percentage_done = activeEmitters[0]['emitter'].percentage_done(steps)
            print str(percentage_done) + '% is emitted, ' + str(len(p)) + ' electrons in the system'
    
        m.one_step()
        steps = steps + 1
        # breaks loop
        if setup.steps != -1 and setup.steps < steps:
            moreIterations = 0            	
    data = [timeLastDeposit, steps, posAll, allAbsorbed, allEmitted]
    return data	
	
	
	
#----=======#######RUN SIMULATION#######=======----	
dataAll.append(emit(data))

file = open("Data.pickle","w")
pickle.dump(dataAll, file)
file.close()

file = open("Absorbed.pickle","w")
pickle.dump(allAbsorbed['nr'], file)
file.close()

#----=======#######DATA PLOTS#######=======----
allAbsorbed = data[3]
allEmitted = data[4]

#pulseTemp = analyze_tools.smoothAbsorptionTimeProfile(dataAll[-1], 100)
#pylab.clf()
#pylab.plot(pulseTemp)
#pylab.ylim([0,.12])
#pylab.savefig("Pulse.png")
	
pylab.clf()
plot_tools.plotSampling(nr)
pylab.ylabel("#")
pylab.xlabel("#")
pylab.savefig("Eletrons_In_Gap_vs_Iter.png")

if 1: #absorber stuff
    pylab.clf()
    plot_tools.plotSampling(allAbsorbed['nr'])
    pylab.ylabel("#")
    pylab.xlabel("#")
    pylab.savefig("Absorbed_vs_Iter.png")

    pylab.clf()
    p = analyze_tools.determineAbsorptionTimeProfile(allAbsorbed['nr'],100)
    pylab.ylabel("#")
    pylab.xlabel("t")
    pylab.plot(p)
    pylab.savefig("Absorbed_vs_Time.png")

#    current = analyze_tools.determineCurrentDensity(allAbsorbed['nr'])
#    pylab.clf()
#    plot_tools.plotSampling(current)
#    pylab.ylabel("I")
#    pylab.xlabel("#")
#    pylab.savefig("Absorbed_Current_vs_Iter.png")
    
    absorptionPlane = analyze_tools.determineAbsorptionProfilePlane(allAbsorbed['coordinates'], 100)    
    pylab.clf()
    pylab.axis('equal')
    pylab.pcolormesh(absorptionPlane)
    pylab.xticks( [0,50,100], (str(-setup.x/constants.nm), '0', str(setup.x/constants.nm)))
    pylab.ylabel("nm")
    pylab.yticks( [0,50,100], (str(-setup.y/constants.nm), '0', str(setup.y/constants.nm)))
    pylab.xlabel("nm")
    pylab.colorbar()
    pylab.savefig("Absorbed_In_Plane.png")

    absorptionLine = analyze_tools.determineAbsorptionProfileLine(allAbsorbed['coordinates'], 50)
    pylab.clf()
    x = scipy.array(range(0, len(absorptionLine)))
    x = x*setup.x/(len(absorptionLine)*constants.nm)
    pylab.plot(x, absorptionLine)
    pylab.xlabel("nm")
    pylab.savefig("Absorbed_Radial_Profile.png")

if 1: # emitter stuff
    pylab.clf()
    plot_tools.plotSampling(allEmitted['nr'])
    pylab.ylabel("#")
    pylab.xlabel("#")
    pylab.savefig("Emitted_vs_Iter.png")

    pylab.clf()
    p = analyze_tools.determineAbsorptionTimeProfile(allEmitted['nr'],100)
    pylab.ylabel("#")
    pylab.xlabel("t")
    pylab.plot(p)
    pylab.savefig("Emitted_vs_Time.png")

#    current = analyze_tools.determineCurrentDensity(allEmitted['nr'])
#    pylab.clf()
#    plot_tools.plotSampling(current)
#    pylab.ylabel("I")
#    pylab.xlabel("#")
#    pylab.savefig("Emitted_Current_vs_Iter.png")

    emissionPlane = analyze_tools.determineAbsorptionProfilePlane(allEmitted['coordinates'], 100)    
    pylab.clf()
    pylab.axis('equal')
    pylab.pcolormesh(emissionPlane)
    pylab.xticks( [0,50,100], (str(-setup.x/constants.nm), '0', str(setup.x/constants.nm)))
    pylab.ylabel("nm")
    pylab.yticks( [0,50,100], (str(-setup.y/constants.nm), '0', str(setup.y/constants.nm)))
    pylab.xlabel("nm")
    pylab.colorbar()
    pylab.savefig("Emitted_In_Plane.png")

    emissionLine = analyze_tools.determineAbsorptionProfileLine(allEmitted['coordinates'], 50)
    pylab.clf()
    x = scipy.array(range(0, len(emissionLine)))
    x = x*setup.x/(len(emissionLine)*constants.nm)
    pylab.plot(x, emissionLine)
    pylab.xlabel("nm")
    pylab.savefig("Emitted_Radial_Profile.png")

