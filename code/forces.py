import force_calculators
reload(force_calculators)

#import scipy

class Forces:
    def __init__(self, N, on=1):
#		self._externalPotential = 0
        self._forces = {}
        self._forces["InterParticle"] = force_calculators.Coulomb(N, on)
        return
    
    def set_constant_external_field(self, N, E):
        #returns true if any of the elements are defined
        if not E.any():
            self._forces.pop("ConstantExternal")
        else:
            self._forces["ConstantExternal"] = force_calculators.ExternalConstantField(E, N)
        return

    def set_linear_external_field(self, N, E, radial=0):
        #returns true if any of the elements are defined
        if not E.any():
            self._forces.pop("LinearExternal")
        else:
            self._forces["LinearExternal"] = force_calculators.ExternalLinearField(E, N, radial)
        return
        
        
    def set_quadratic_external_field(self, N, E, radial=0):
        #returns true if any of the elements are defined
        if not E.any():
            self._forces.pop("QuadraticExternal")
        else:
            self._forces["QuadraticExternal"] = force_calculators.ExternalQuadraticField(E, N, radial)
        return


                        
    def get_F(self):
        keys = self._forces.keys()
        F = self._forces[keys[0]].get_F()
        for key in keys[1:]:			
            F = F+self._forces[key].get_F()
        return F
        
    def get_single_F(self):
        keys = self._forces.keys()
        F = self._forces[keys[0]].get_single_F()
        for key in keys[1:]:
            F = F+self._forces[key].get_single_F()
        return F
        
    def remove(self, i):
        for key in self._forces.keys():
            self._forces[key].remove(i)
        return

    def add(self):
        for key in self._forces.keys():
            self._forces[key].add()
        return

    def determine_forces(self, R):
        for key in self._forces.keys():
            self._forces[key].determine_forces(R)
        return
        
    def determine_forces_single(self, R, i):
        for key in self._forces.keys():
            self._forces[key].determine_forces_single(R, i)
        return

    def set_box(self, box):
        for key in self._forces.keys():
            self._forces[key].set_box(box)
        return
