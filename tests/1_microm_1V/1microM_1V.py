import scipy
from scipy import stats
import pylab
import random
import pickle

import particles
reload(particles)

import movers as mover
reload(mover)

import emitters as emitter
reload(emitter)

import boundaries as boundary
reload(boundary)

import helper_functions
reload(helper_functions)

import constants
reload(constants)

import plot_tools
reload(plot_tools)

import analyze_tools
reload(analyze_tools)

import Input
reload(Input)



def add_iterations(iterations, data):
	timeLastDeposit = data[0]
	steps = data[1]
	posAll = data[2]
	phase = data[3]
	absorbed = data[4]
	emitted = data[5]
	
	for i in xrange(0,iterations):
		# impose boundary conditions
		for boundary in boundaries:
			absorbedStepAll = boundary.checkAll(p)
			# keep track on how many atoms that have been absorbed
			if absorbedStepAll:
				for absorbedStep in absorbedStepAll:
					absorbed['cordinates'].append(absorbedStep)
		absorbed['nr'].append(len(absorbed['cordinates']))		

		totalTime.append(totalTime[-1] + Input.timeStep)
		
		if i%Input.sampling == 0:
			print "Percent completed: "+str(int(i/float(iterations)*100))+'%'
			nr.append(len(p))
			if posAll is not None:
				pos = p.getR()
				posAll.append(pos)

		# emit new particles
		nrParticlesBefore = len(p)
		if Input.timeDeposite < timeLastDeposit:
			for emitter in emitters:
				emitter.addParticles(p, 1)
			timeLastDeposit = 0.

		emitted['nr'].append(len(p)-nrParticlesBefore)
		
		# stop if the gab is not containing any particles
#		if 	len(p) == 0:
#			break
		# stop when X particles have been absorbed
#		if 20000 < absorbed['nr'][-1]:
#			break
		# stop if any particle is moving backwards (Should be checked if 2D simulation)
#		if i and (min(v[:,2]) < 0.):
#			break

		# evolve the system one discrete time step
		m.oneStep()
		steps = steps + 1
		timeLastDeposit = timeLastDeposit + Input.timeStep
	
		if iterations-5000 < i and phase is not None:
			pos_i = p.getR()
			vel_i = p.getV()
			acc_i = p.getA()
			for j in xrange(len(pos_i)):
				phase['r'].append(pos_i[j,-1])
				phase['v'].append(vel_i[j,-1])
				phase['a'].append(acc_i[j,-1])
			
	data = [timeLastDeposit, steps, posAll, phase, absorbed, emitted]
	return data


	
#----=======#######SINGLE PULSE#######=======----	
def single_pulse(data):
	timeLastDeposit = data[0]
	steps = data[1]
	posAll = data[2]
	absorbed = data[3]
	emitted = data[4]
		
	# gives a pulse length t containing x particles
#	emissionEvents = helper_functions.make_gaussian_distribution(1.e-11,50)	
#	emissionEvents = helper_functions.make_square_distribution(1.e-11,1500)

#	emissionEvents = helper_functions.make_square_distribution(3.e-11,3000)
	emissionEvents = helper_functions.make_square_distribution(3.e-13,3000)
	
	iEmission = 1
#	extraSteps = 100
	extraSteps = 40000

	while 1:
		# impose boundary conditions
		for boundary in boundaries:
			absorbedStepAll = boundary.check_all(p)
			# keep track on how many atoms that have been absorbed
			if absorbedStepAll:
				for absorbedStep in absorbedStepAll:
					absorbed['cordinates'].append(absorbedStep)
		absorbed['nr'].append(len(absorbed['cordinates']))

		totalTime.append(totalTime[-1] + Input.timeStep)

		if steps%Input.sampling == 0:
			nr.append(len(p))
			if posAll is not None:
				pos = p.get_R()
				posAll.append(pos)

		# emit new particles
		nrParticlesBefore = len(p)
		try:
			if emissionEvents[iEmission] < timeLastDeposit:
				for emitter in emitters:
					emitter.add_particles(p, 1)
				timeLastDeposit = 0.
				iEmission = iEmission + 1
#				print iEmission
		# emission process is over
		except IndexError:
			if 	extraSteps < 0:
				break
			else:
				extraSteps = extraSteps-1

		emitted['nr'].append(len(p)-nrParticlesBefore)

		m.one_step()
		steps = steps + 1
		timeLastDeposit = timeLastDeposit + Input.timeStep
	
	data = [timeLastDeposit, steps, posAll, absorbed, emitted]
	return data	




#----=======#######DRIVING CODE#######=======----

dataAll = []
pulsesAll = []
#for kk in scipy.array([3.0, 4.0, 5.0]):
for kk in scipy.array([0.0]):
#	for jj in scipy.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ]):
	for jj in scipy.array([0.5]):
#		for ii in kk+jj+scipy.array([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09 ]):
		for ii in kk+jj+scipy.array([0.0 ]):
			print ii
			
#			r = scipy.zeros((Input.particles,Input.dim))
#			p = Particles.Particles(len(r),Input.dim, Input.coulombPotetial)
			p = particles.Particles(0)

			#----=======#######EXTERNAL POTENTIAL#######=======----
			E_z = ii/(2.*Input.z)
			CosntantE = scipy.array([0., 0., E_z])
			
			if Input.CosntantE.any():
		#		p.setExternalPotential(Input.CosntantE, "Constant")
				p.set_external_potential(CosntantE, "Constant")

			if Input.LinearE.any():
				p.set_external_potential(Input.LinearE, "Linear")

			if Input.RadialLinearE:
				p.set_external_potential(scipy.array([Input.RadialLinearE, 0., 0.]), "RadialLinear")

			emitters = []
			boundaries = []
#			#----=======#######2D#######=======----
#			if Input.dim == 2:
#				p.setBox(scipy.array([Input.x, Input.y])*Input.periodic)
#			#----=======#######EMITTER#######=======----
#				emitters.append(Emitter.Line(-Input.y))
#			#----=======#######BOUNDARIES#######=======----		
#				boundaries.append(Boundary.Min('Reflecting', 'x', -Input.x))
#				boundaries.append(Boundary.Max('Reflecting', 'x',  Input.x))
#				boundaries.append(Boundary.Min('Reflecting', 'y', -Input.y))
#				boundaries.append(Boundary.Max('Absorbing',  'y',  Input.y))
#
#			#----=======#######3D#######=======----
			if 1:
				p.set_box(scipy.array([Input.x, Input.y, Input.z])*Input.periodic)
			#----=======#######EMITTER#######=======----	
				if Input.shapeDeposit == 'rectangel':
					emitters.append(Emitter.Plane(-Input.zDeposit, Input.xDeposit, -Input.xDeposit, Input.yDeposit, -Input.yDeposit))
				elif Input.shapeDeposit == 'pattern':
					print "pattern"
			##		emitters.append(Emitter.Plane(-Input.zDeposit,  Input.xDeposit,     Input.xDeposit*1./3.,  Input.yDeposit,     Input.yDeposit*1./3.))
			##		emitters.append(Emitter.Plane(-Input.zDeposit, -Input.xDeposit*1./3., -Input.xDeposit,     Input.yDeposit,     Input.yDeposit*1./3.))
			##		emitters.append(Emitter.Plane(-Input.zDeposit, -Input.xDeposit*1./3., -Input.xDeposit,    -Input.yDeposit*1./3., -Input.yDeposit   ))
			##		emitters.append(Emitter.Plane(-Input.zDeposit,  Input.xDeposit,     Input.xDeposit*1./3., -Input.yDeposit*1./3., -Input.yDeposit   ))
					emitters.append(emitter.Plane(-Input.zDeposit,  Input.xDeposit,     Input.xDeposit/3.,  Input.yDeposit/3., -Input.yDeposit/3.))
					emitters.append(emitter.Plane(-Input.zDeposit, -Input.xDeposit/3., -Input.xDeposit,     Input.yDeposit/3., -Input.yDeposit/3.))
					
					emitters.append(emitter.Plane(-Input.zDeposit,  Input.xDeposit/3., -Input.xDeposit/3., -Input.yDeposit/3., -Input.yDeposit   ))
					emitters.append(emitter.Plane(-Input.zDeposit,  Input.xDeposit/3., -Input.xDeposit/3.,  Input.yDeposit,    Input.yDeposit/3.))
					
			##		emitters.append(Emitter.Plane(-Input.zDeposit,  Input.xDeposit/3., -Input.xDeposit/3.,  Input.yDeposit/3., -Input.yDeposit/3.))
				elif Input.shapeDeposit == 'circle':
			#		emitters.append(Emitter.PlaneSpaceCharge(-Input.zDeposit, Input.rDeposit))
					emitters.append(emitter.PlaneRoundProfile(-Input.zDeposit, Input.rDeposit, 'Circle'))

				elif Input.shapeDeposit == 'repeat':
					file = open("emission.pickle", "r")
					completePuls = pickle.load(file)
					file.close()
					emitters.append(emitter.RepeatEmission(completePuls))
					
			#----=======#######BOUNDARIES#######=======----
				boundaries.append(boundary.Min('Reflecting', 'x', -Input.x))
				boundaries.append(boundary.Max('Reflecting', 'x',  Input.x))
				boundaries.append(boundary.Min('Reflecting', 'y', -Input.y))
				boundaries.append(boundary.Max('Reflecting', 'y',  Input.y))
				boundaries.append(boundary.Min('Reflecting', 'z', -Input.z))
				boundaries.append(boundary.Max('Absorbing',  'z',  Input.z))

			#	boundaries.append(boundary.PointInPlane('Reflecting', 'z', 100.*Input.nm, -Input.z/2.))

				
			#----=======#######DATA STORAGE#######=======----
			# to ensure that a particle is initially deposited
			timeLastDeposit = Input.timeLastDeposit
			steps = 0

			posAll = []

			phase = None
			if phase is not None:
				phase['r'] = []
				phase['v'] = []
				phase['a'] = []

			absorbed = {}
			absorbed['cordinates'] = []
			absorbed['nr'] = []

			emitted = {}
			emitted['nr'] = []

			nr = []
			totalTime = list([0.0])

			m = mover.VelocityVerlet(p, Input.timeStep)
			#data = [timeLastDeposit, steps, posAll, phase, absorbed, emitted]
			#data = addIterations(Input.iterations, data)
			data = [timeLastDeposit, steps, posAll, absorbed, emitted]
			dataAll.append(single_pulse(data))

			pulseTemp = analyze_tools.smoothAbsorptionTimeProfile(dataAll[-1], 100)
			pulsesAll.append(pulseTemp)
			pylab.clf()
			pylab.plot(pulseTemp)
			pylab.ylim([0,.12])
			pylab.savefig("pulse"+str(ii)+".eps")
			pylab.savefig("pulse"+str(ii)+".png")
	

#----=======#######DATA ANALYSIS#######=======----
absorbed = data[3]
current = analyze_tools.determineCurrentDensity(absorbed['nr'])
absorptionPlane = analyze_tools.determineAbsorptionProfilePlane(absorbed['cordinates'], 100)
absorptionLine = analyze_tools.determineAbsorptionProfileLine(absorbed['cordinates'], 50)

pylab.clf()
pylab.axis('equal')
pylab.pcolormesh(absorptionPlane)
pylab.xticks( [0,50,100], (str(-Input.x/constants.nm), '0', str(Input.x/constants.nm)))
pylab.ylabel("nm")
pylab.yticks( [0,50,100], (str(-Input.y/constants.nm), '0', str(Input.y/constants.nm)))
pylab.xlabel("nm")
pylab.colorbar()
pylab.savefig(str(Input.V_z_constant)+"V_z_constant_"+str(2.*Input.z/constants.nm)+"nm_Absorbed_Plane_1.png")
#pylab.savefig('1.png')
pylab.clf()
x = scipy.array(range(0, len(absorptionLine)))
x = x*Input.x/(len(absorptionLine)*constants.nm)
pylab.plot(x, absorptionLine)
pylab.xlabel("nm")
pylab.savefig(str(Input.V_z_constant)+"V_z_constant_"+str(2.*Input.z/constants.nm)+"nm_Absorbed_Radial_1.png")
#pylab.savefig('2.png')

pylab.clf()
plot_tools.plotSampling(nr)
pylab.ylabel("#")
pylab.savefig(str(Input.V_z_constant)+"V_z_constant_"+str(2.*Input.z/constants.nm)+"nm_Nr_vs_Iter_1.png")

pylab.clf()
plot_tools.plotSampling(current)
pylab.ylabel("I")
pylab.savefig(str(Input.V_z_constant)+"V_z_constant_"+str(2.*Input.z/constants.nm)+"nm_Current_vs_Iter_1.png")

pylab.clf()
plot_tools.plotSampling(absorbed['nr'])
pylab.ylabel("#")
pylab.savefig(str(Input.V_z_constant)+"V_z_constant_"+str(2.*Input.z/constants.nm)+"nm_"+str(Input.rDeposit/constants.nm)+"nm_Absorbed_vs_Iter_1.png")

file = open("Absorbed_Emitter_"+str(Input.rDeposit/constants.nm)+"nm.pickle","w")
pickle.dump(absorbed['nr'], file)
file.close()


pylab.clf()
p = analyze_tools.determineAbsorptionTimeProfile(absorbed['nr'],50)
pylab.ylabel("#")
pylab.ylabel("t")
pylab.plot(p)
pylab.savefig(str(Input.V_z_constant)+"V_z_constant_"+str(2.*Input.z/constants.nm)+"nm_"+str(Input.rDeposit/constants.nm)+"nm_Absorbed_vs_Time.png")


