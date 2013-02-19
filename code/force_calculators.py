import numpy

import potentials.coulomb
reload(potentials.coulomb)

class ForceCalculator:
    def __init__(self, N):
        self._N = N
        self._F = numpy.zeros((self._N, 3))
        self._singleF = numpy.zeros(3)
        return
    
    def get_F(self):
        return self._F

    def get_single_F(self):
        return self._singleF
        
    def remove(self, i):
        self._N = self._N-1
        self._F = numpy.reshape(numpy.append(self._F[:i], self._F[i+1:]),
                                             (self._N, 3))
        return

    def add(self):
        self._N = self._N+1
        self._F = numpy.reshape(numpy.append(self._F, numpy.zeros(3)),
                                             (self._N, 3))
        return

    def determine_forces(self):
        pass
        return

    def determine_forces_single(self):
        pass
        return
        
    def set_box(self, box):
        self.box_ = box
        return

#------=====######EXTERNAL######=====------

class ExternalConstantField(ForceCalculator):
    def __init__(self, E0, N):
        ForceCalculator.__init__(self, N)
        self._F_E = E0*1.6021e-19		 #E*e in si units
        return
                
    def determine_forces(self, R):
        assert len(R) == self._N
        tempF = self._singleF

        for i in xrange(0, self._N):
            self.determine_forces_single(R, i)
            self._F[i] = self._singleF
        self._singleF = tempF
        return
        
    def determine_forces_single(self, R, i):
        self._singleF = self._F_E
        return

class ExternalLinearField(ForceCalculator):
    def __init__(self, E0, N, radial=None):
        ForceCalculator.__init__(self, N)
        self.radial_ = radial
        self._F_E = E0*1.6021e-19		 #E*e in si units
        return
                
    def determine_forces(self, R):
        assert len(R) == self._N
        tempF = self._singleF

        for i in xrange(0, self._N):
            self.determine_forces_single(R, i)
            self._F[i] = self._singleF
        self._singleF = tempF
        return

    def determine_forces_single(self, R, i):
        if self.radial_:
            r2 = R[i][0]**2 + R[i][1]**2
            r = numpy.sqrt(r2)
            f = r*self._F_E[0]
            #handles singularity at 1/0	
            if r == 0.0:
                self._singleF[0] = 0.0
                self._singleF[1] = 0.0
            else:
                self._singleF[0] = -R[i][0]/r * f 
                self._singleF[1] = -R[i][1]/r * f 
        else:
            #x
            if self._F_E[0]:
                self._singleF[0] = -R[i][0]*self._F_E[0] 
            else:
                self._singleF[0] = 0.0
            #y
            if self._F_E[1]:
                self._singleF[1] = -R[i][1]*self._F_E[1] 
            else:
                self._singleF[1] = 0.0
            #z
            if self._F_E[2]:
                self._singleF[2] = -R[i][2]*self._F_E[2] 
            else:
                self._singleF[2] = 0.0
        return
        
class ExternalQuadraticField(ForceCalculator):
    def __init__(self, E0, N, radial=None):
        ForceCalculator.__init__(self, N, 3)
        self.radial_ = radial
        self._F_E = E0*1.6021e-19		 #E*e in si units
        return
                
    def determine_forces(self, R):
        assert len(R) == self._N
        tempF = self._singleF

        for i in xrange(0, self._N):
            self.determineForcesSingle(R, i)
            self._F[i] = self._singleF
        self._singleF = tempF
        return

    def determine_forces_single(self, R, i):
        if self.radial_:
            r2 = R[i][0]**2 + R[i][1]**2
            r = numpy.sqrt(r2)
            f = r2*self._F_E[0]
            #handles singularity at 1/0	
            if r == 0.0:
                self._singleF[0] = 0.0 
                self._singleF[1] = 0.0 			
            else:
                self._singleF[0] = -R[i][0]/r * f 
                self._singleF[1] = -R[i][1]/r * f 
        else:
            #x
            if self._F_E[0]:
                self._singleF[0] = -sign(R[i][0]) * R[i][0]**2*self._F_E[0] 
            else:
                self._singleF[0] = 0.0
            #y
            if self._F_E[1]:
                self._singleF[1] = -sign(R[i][1]) * R[i][1]**2*self._F_E[1] 
            else:
                self._singleF[1] = 0.0
            #z
            if self._F_E[2]:
                self._singleF[2] = -sign(R[i][2]) * R[i][2]**2*self._F_E[2] 
            else:
                self._singleF[2] = 0.0
        return

#------=====######INTERNAL######=====------

class Coulomb(ForceCalculator):
    def __init__(self, N, on=0):
        ForceCalculator.__init__(self, N)
        pos = numpy.zeros((self._N,3))
        self._potential = potentials.coulomb.Coulomb(pos)
        self._potential.set_prefator(2.30686e-28) #in si units, 2.30686e-28
                                                 #e^2*(4*pi*eps_0)
                                                 #(-1.6021e-19)**2*8.9876e9
        self._potential.set_R0(0.0)				 #to prevent singularity a r=0
        self._on = on
        return

    def set_box(self, box):
        ForceCalculator.set_box(self, box)
        self._potential.set_box(box[0], box[1], box[2])
        return		
        
    def add(self):
        ForceCalculator.add(self)
        self._potential.change_nr_atoms(self._N)
        return
        
    def remove(self, i):
        ForceCalculator.remove(self, i)
        self._potential.change_nr_atoms(self._N)
        return
        
    def determine_forces(self, R):
        assert len(R) == self._N
        if self._on:
            self._potential.set_positions(R)
            self._F = self._potential.get_forces()
        else:
            self._F = numpy.array(R)*0.
        return
        
    def determine_forces_single(self, R, i):
        assert len(R) == self._N
        if self._on:
            self._potential.set_positions(R)
            self._singleF = self._potential.get_forces_single(i)
        else:
            self._singleF = numpy.array(R[i])*0.
        return
