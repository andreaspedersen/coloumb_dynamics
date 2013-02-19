import numpy
from c_potentials import Coulomb_C

class Coulomb:
    def __init__(self, pos=numpy.array([[.0, .0, .0]])):
        self._coulomb = Coulomb_C()
        self._coulomb.init(len(pos))							

        self.set_positions(pos)
        
        self._F = numpy.array(pos)*0.
        return

    def set_prefator(self, prefactor):
        self._coulomb.setPrefactor(prefactor)
        return

    def set_R0(self, R0):
        self._coulomb.setR0(R0)
        return

    def set_box(self, x, y, z):
        self._coulomb.setBox(x, y, z)
        return

    def set_positions(self, pos):
        for i in xrange(0,len(pos)):
            self._coulomb.setPosition(i, float(pos[i,0]), float(pos[i,1]), float(pos[i,2]))
        return

    def get_forces(self):
        self._F = self._F*0.
        for i in xrange(0,len(self._F)):
            self._F[i] = [self._coulomb.getForce(i,0), 
                          self._coulomb.getForce(i,1),
                          self._coulomb.getForce(i,2)]
        return self._F

    def get_forces_single(self, i):
        singleF = numpy.array([self._coulomb.getForceSingle(i,0), 
                               self._coulomb.getForceSingle(i,1),
                               self._coulomb.getForceSingle(i,2)])
        return singleF
        
    def change_nr_atoms(self, nrAtoms):
        self._coulomb.setNrAtoms(nrAtoms)
        self._F = numpy.zeros((nrAtoms,3))
        return		

