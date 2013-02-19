import numpy

import forces
reload(forces)

import constants
reload(constants)

class Particles:
    def __init__(self, N, on=1):
        self._updateForce = 1
        self._updateForceSingle = 1

        self._N = N
        self._coulomb_on = on

        self._box = numpy.zeros(3)
        self._R = numpy.zeros((self._N,3))
        self._V = numpy.zeros((self._N,3))
        self._F = numpy.zeros((self._N,3))
        self._forcesCalc = forces.Forces(self._N, on)

        self._externalE = {}
        self.timeUsed = 0.
        return

    def __len__(self): 
        return self._N  

    def __copy__(self):
        newParticles = Particles(self._N, self._coulomb_on)
        newParticles.set_box(self._box)
        for key in self._externalE.keys():
            newParticles.set_external_potential(self._externalE[key], key)
        newParticles.set_R(self._R)
        newParticles.set_V(self._V)
        return newParticles	
            
    def set_external_potential(self, E, type):
        self._updateForce = 1
        self._updateForceSingle = 1

        self._externalE[type] = E
        if type == 'Constant':
            self._forcesCalc.set_constant_external_field(self._N, self._externalE[type])
        elif type == 'Linear':
            self._forcesCalc.set_linear_external_field(self._N, self._externalE[type])
        elif type == 'RadialLinear':
            radial = 1
            self._forcesCalc.set_linear_external_field(self._N, self._externalE[type], radial)
        elif type == 'Quadratic':
            self._forcesCalc.set_quadratic_external_field(self._N, self._externalE[type])
        elif type == 'RadialQuadratic':
            radial = 1
            self._forcesCalc.set_quadratic_external_field(self._N, self._externalE[type], radial)
        else:
            print "STOP: Not a known external potential: " + str(type) 
            stop
        return
        
    def set_R(self, R):
        self._updateForce = 1
        self._updateForceSingle = 1

        assert self._N == len(R)
        self._R = numpy.array(numpy.reshape(R,(self._N,3)))
        return
    
    def set_V(self, V):
        self._updateForce = 1
        self._updateForceSingle = 1

        assert self._N == len(V)
        self._V = numpy.array(numpy.reshape(V,(self._N,3)))
        return

    def set_box(self, box):
        self._updateForce = 1
        self._updateForceSingle = 1

        assert 3 == len(box)
        self._box = box
        self._forcesCalc.set_box(self._box)
        return
    
    def get_box(self):
        return self._box
    
    def get_R(self):
        return self._R

    def get_V(self):
        return self._V

    def get_F(self):
        if self._updateForce == 1:
            self._forcesCalc.determine_forces(self._R)
            self._F = self._forcesCalc.get_F()
            self._updateForce = 0
        return self._F
    
    def get_single_F(self, i):
        if self._updateForceSingle == 1:
            self._forcesCalc.determine_forces_single(self._R, i)
            self.single_F = self._forcesCalc.get_single_F()
            self._updateForceSingle = 0
        return self.single_F

    def get_A(self):
        return self.get_F()/constants.me
    
    def remove(self, i):
        self._updateForce = 1
        self._updateForceSingle = 1

        self._N = self._N-1
        # last element
        if i == self._N:
            self._R = numpy.reshape(self._R[:i],(self._N,3))
            self._V = numpy.reshape(self._V[:i],(self._N,3))
            self._F = numpy.reshape(self._F[:i],(self._N,3))
        else:
            self._R = numpy.reshape(numpy.append(self._R[:i], self._R[i+1:]),(self._N,3))
            self._V = numpy.reshape(numpy.append(self._V[:i], self._V[i+1:]),(self._N,3))
            self._F = numpy.reshape(numpy.append(self._F[:i], self._F[i+1:]),(self._N,3))

        self._forcesCalc.remove(i)
        return

    def add(self, r, v=None, f=None):
        self._updateForce = 1
        self._updateForceSingle = 1

        if v is None:
            v = numpy.zeros(3)
        if f is None:
            f= numpy.zeros(3)
            
        self._N = self._N+1
        self._R = numpy.reshape(numpy.append(self._R, r),(self._N,3))
        self._V = numpy.reshape(numpy.append(self._V, v),(self._N,3))
        self._F = numpy.reshape(numpy.append(self._F, f),(self._N,3))
        self._forcesCalc.add()
        return
        
