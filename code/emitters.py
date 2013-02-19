import numpy
import random
import copy
import pickle

import setup
reload(setup)

import helper_functions
reload(helper_functions)


class Emitter:
    def __init__(self, zMin=None, xMax=None, xMin=None, yMax=None, yMin=None, filename="Emitted.pickle"):

        if setup.emissionShape == 'circle':
            self._rMax = xMax
            self._rMax2 = self._rMax**2
            self._z = zMin

            if xMin == None:
                self._x = 0.
            else:
                self._x = xMin   
            if yMax == None:
                self._y = 0.
            else:
                self._y = yMax

        elif setup.emissionShape == 'rectangel':
            self._xMin = xMin
            self._xMax = xMax
            self._dx = xMax-xMin
            self._yMin = yMin
            self._yMax = yMax
            self._dy = yMax-yMin
            self._rMax2 = None
            self._z = zMin

        self._completePuls = []        
        self._filename = filename
        random.seed()
        return
        
    def __del__(self):
        self.save_complete_pulse()
        
    def attempt_add_particle(self, particlesTest):
        add_pos = None
        add_vel = None
        emit = 0 
        
        z = self._z
        #rectangle
        if self._rMax2 is None:
            x = self._xMin+random.random()*self._dx
            y = self._yMin+random.random()*self._dy
        #circle
        else:
            r2 = self._rMax2+1.
            while(self._rMax2 < r2):
                x = numpy.sign(.5-random.random())*random.random()*self._rMax
                y = numpy.sign(.5-random.random())*random.random()*self._rMax
                r2 = x**2+y**2
            x = x + self._x
            y = y + self._y

        posBelow = numpy.array([x, y, z-1.e-9])
        pos = numpy.array([x, y, z])
               
        # could be important when determining if a particle with a finite vel
        # can escape the field at the emitter 
        vx = 0.
        vy = 0.
        vz = 0.
        vel = numpy.array([vx, vy, vz])

        particlesTest.add(posBelow)
        #cannot use -1 as index for the last element!
        indexLast = len(particlesTest)-1
        forceBelowPlane = particlesTest.get_single_F(indexLast)
        particlesTest.remove(len(particlesTest)-1)

        if 0.0 < forceBelowPlane[-1]:
            particlesTest.add(pos)
            #cannot use -1 as index for the last element!
            indexLast = len(particlesTest)-1
            forceInPlane = particlesTest.get_single_F(indexLast)
            #no particle should be deposited, ensures that force always will be negative
            if numpy.isnan(forceInPlane[-1]) == True:
                forceZ = -1.
            else:
                forceZ = min(forceInPlane[-1], forceBelowPlane[-1])

            #determine if a particle should be added
            if (0.0 < forceZ):
                add_vel = numpy.array([vx, vy, vz])
                add_pos = pos
                self._completePuls.append([add_pos, add_vel])

        return add_pos, add_vel

    def add_particle(self, particlesTest):
        emitted = None
        tries = setup.emissionNrAttempts
        while(tries):
            add_pos, add_vel = self.attempt_add_particle(particlesTest)
            if add_pos != None:
                emitted = [add_pos, add_vel]
                break
            else:
                tries = tries-1
        return emitted
    
    def get_complete_pulse(self):
        return self._completePuls

    def save_complete_pulse(self):
        file = open(self._filename,"w")
        pickle.dump(self._completePuls, file)
        file.close()
        return
        
    def percentage_done(self, steps):
        return

class PlaneSpaceCharge(Emitter):
    def __init__(self, zMin, xMax, xMin=None, yMax=None, yMin=None):
        Emitter.__init__(self, zMin, xMax, xMin, yMax, yMin)
        return

    def add_particles(self, particles):
        allEmitted = []
        particlesTest = copy.copy(particles)
        while(1):
            emitted = self.add_particle(particlesTest)
            if emitted == None:
                break
            else:
                allEmitted.append(emitted)

        # nothing was emitted
        if len(allEmitted) == 0:
            allEmitted = 0
        return allEmitted

    def percentage_done(self, steps):
        return int(steps / (setup.steps/100.))

    
#    def add_particles(self, particles):
#        emitted = []
#        particlesTest = copy.copy(particles)
#        tries = setup.emissionNrAttempts
#        while(tries):
#            add_pos, add_vel = self.attempt_add_particle(particlesTest)
#            if add_pos != None:
#                particles.add(add_pos, add_vel)
#                tries = setup.emissionNrAttempts
#                emitted.append([add_pos, add_vel])
#            else:
#                tries = tries-1
#        # nothing was emitted
#        if len(emitted) == 0:
#            emitted = 0  
#        return emitted


class PlanePulse(Emitter):
    def __init__(self, zMin, xMax, xMin=None, yMax=None, yMin=None, eTime=None, eNr=None, eShape=None):
        Emitter.__init__(self, zMin, xMax, xMin, yMax, yMin)

        if eTime == None:
            eTime = setup.emissionDuration
        if eNr == None:
            eNr = setup.emissionNrElectrons
        if eShape == None:
            ePulseShape = setup.emissionPulseShape

        self._emissionLastItertion = 0
        self._iEmission = 0
        self._timePassedSinceDeposit = 0

        if ePulseShape == 'gaussian':
            self._emissionEvents = helper_functions.make_gaussian_distribution(eTime, eNr)
        elif  ePulseShape == 'square':
            self._emissionEvents = helper_functions.make_square_distribution(eTime, eNr)
        else:
            print "The emission pulse shape (emissionPulseShape) is not avaliable"
        return
    

    def add_particles(self, particles, timeStep=None):

        if timeStep == None:
            timeStep = setup.timeStep
        self._timePassedSinceDeposit = self._timePassedSinceDeposit + timeStep

        if self._iEmission < len(self._emissionEvents):
            timeNextDeposit = 0 #self._emissionEvents[self._iEmission]
            allEmitted = []
            blocked = 0
            particlesTest = copy.copy(particles)

            while ((timeNextDeposit < self._timePassedSinceDeposit) and 
                   (self._iEmission < len(self._emissionEvents))):
                timeNextDeposit = self._emissionEvents[self._iEmission]            
                if not blocked:
                    emitted = self.add_particle(particlesTest)
                    if emitted == None:
                        blocked = 1
                    else:
                        allEmitted.append(emitted)
                    self._timePassedSinceDeposit = self._timePassedSinceDeposit - timeNextDeposit
                    self._iEmission = self._iEmission + 1
        else:
            allEmitted = -1
        return allEmitted

    def percentage_done(self, steps):
        return int(self._iEmission / (len(self._emissionEvents)/100.))
            

class Plane(Emitter):
    def __init__(self, zMin, xMax, xMin=None, yMax=None, yMin=None):
        Emitter.__init__(self, zMin, xMax, xMin, yMax, yMin)

    def add_particles(self, particles, addNrParticles):
        emitted = []
        particlesTest = copy.copy(particles)
        while(addNrParticles):
            tries = setup.emissionNrAttempts
            while(tries):
                add_pos, add_vel = self.attempt_add_particle(particlesTest)
                if add_pos != None:
                    particles.add(add_pos, add_vel)
                    addNrParticles = addNrParticles - 1
                    emitted.append([add_pos, add_vel])
                    break
                else:
                    tries = tries-1
            if tries == 0:
                addNrParticles = 0
#                print "give up!"
            # nothing was emitted
            if len(emitted) == 0:
                emitted = 0
        return emitted
        
class RepeatEmission(Emitter):
    def __init__(self, filename=None):
        Emitter.__init__(self)
        if filename == None:
            filename = self._filename
        file = open(filename, "r")
        self._completePuls = pickle.load(file)
        file.close()

    def add_particles(self, particles, addNrParticles):
        particlesLeftToAdd = addNrParticles
        while(particlesLeftToAdd):
            nextEmission = self._completePuls.pop(0)
            try:
                len(nextEmission[0])
                particles.add(nextEmission[0],nextEmission[1])
            except TypeError:
                pass
            particlesLeftToAdd = particlesLeftToAdd-1
        return
 
 
 
class LineSpaceCharge:
    def __init__(self, zMin, xMax, xMin):
        self._z = zMin
        self._xMax = xMax
        self._xMin = xMin
        self._dx = xMax-xMin

        return
        
    def add_particles(self, particles, addTries):
        posAll = particles.getR()
        particlesTest = copy.copy(particles)
        tries = addTries
        while(tries):
            z = self._z
            add = 0
            x = self._xMin+random.random()*self._dx

            posBelow = numpy.array([x, z-1.e-9])
            pos = numpy.array([x, z])
                    
            vx = 0.#abs(random.gauss(0,1))*1.e2
            vz = 0.#abs(random.gauss(0,1))*1.e4#abs(random.gauss(0,1))*93776.#0.#velocity giving a E_kin = E_pot(1V)
            vel = numpy.array([vx, vz])

            particlesTest.add(posBelow, vel)
            #cannot use -1 as index for the last element!
            indexLast = len(particlesTest)-1
            forceBelowLine = particlesTest.get_single_F(indexLast)
            particlesTest.remove(len(particlesTest)-1)

            if 0.0 < forceBelowLine[-1]:
                particlesTest.add(pos, vel)
                #cannot use -1 as index for the last element!
                indexLast = len(particlesTest)-1
                forceInLine = particlesTest.get_single_F(indexLast)
                #no particle should be deposited, ensures that forceZ always will be negative
                if numpy.isnan(forceInLine[-1]) == True:
                    forceZ = -1.
                else:
                    #picking out the most negative force
                    if forceInLine[-1] < forceBelowLine[-1]:
                        force = forceInLine
                    else:
                        force = forceBelowLine

                #determine if a particle should be added
                if (0.0 < force[-1]):
                    add = 1
                    particles.add(pos, vel)
            #account for how many unsuccessfull attempts have been carried out in row
            if add == 0:
                tries = tries-1
            else:
                tries = addTries
        return
 
