import numpy

class Mover:
    def __init__(self, particles, dt):
        self._particles = particles
        self._dt = dt
        v = self._particles.get_V()
        a = self._particles.get_A()
        
    def one_step(self):
        pass
        return
    
    def series_of_steps(self, nr):
        for i in xrange(0, nr):
            self.one_step()
        return
        
    def a_nan_in_list(self, list2D):
        aNan = 0
        for a in list2D:
            if numpy.isnan(a[0])  or numpy.isnan(a[1]):# or numpy.isnan(a[2]):
                stop
        return aNan


class VelocityVerlet(Mover):
    def __init__(self, particles, dt):
        Mover.__init__(self, particles, dt)
        return

    def one_step(self):
        r_0 = self._particles.get_R()
        v_0 = self._particles.get_V()
        a_0 = self._particles.get_A()

        self._particles.set_R(r_0+self._dt*v_0+self._dt**2/2.*a_0)
        a_1 = self._particles.get_A()
        self.a_nan_in_list(a_0)
        self.a_nan_in_list(a_1)
        self._particles.set_V(v_0+(a_0+a_1)/2.*self._dt)
    
        return
    
    
