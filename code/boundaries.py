from shapes import Line

import constants
reload(constants)

class Boundary(Line):
    def __init__(self, type, axis, x0, y0, slope):
        self._type = type

        if   axis == 'X' or axis == 'x' or axis == 0:
            self._axis = 0
        elif axis == 'Y' or axis == 'y' or axis == 1:
            self._axis = 1
        elif axis == 'Z' or axis == 'z' or axis == 2:
            self._axis = 2
        else:
            print axis
            "Unkown axis !!"
            stop
        Line.__init__(self, x0, y0, slope)
        return
    
    def check_all(self, particles):
        absorbed = 0
        if self._type == 'Absorbing':
            absorbed = self.absorbing_boundary(particles)
        elif self._type == 'Reflecting':
            self.reflecting_boundary(particles)
        elif self._type == 'Periodic':
            self.periodic_boundary(particles)
        else:
            print "unkown boundary type: "+str(self._type)
        return absorbed
        
    def absorbing_boundary(self, particles):
        pos = particles.get_R()
        absorbed = []
        # has to from (nr atoms-1) to -1
        for i in xrange(len(pos)-1,0-1,-1):
            if self.crossed_boundary(pos[i]):
                absorbed.append(pos[i])
                particles.remove(i)
        if (len(absorbed) == 0):
            absorbed = 0
        return absorbed
    
    def reflecting_boundary(self, particles):
        pos = particles.get_R()
        vel = particles.get_V()
        
        for i in xrange(0, len(pos)):
            if self.crossed_boundary(pos[i]):
                pos[i], vel[i] = self.bounce(pos[i], vel[i])

        particles.set_R(pos)
        particles.set_V(vel)
        return
    
    def periodic_boundary(self, particles):
        pos = particles.getR()
        box = particles.getBox()

        for i in xrange(0,len(pos)):
            if self.crossed_boundary(pos[i]):
                pos[i] = self.flip(pos[i], box)
        particles.setR(pos)
        return

class Min(Boundary):
    def __init__(self, type, axis, x0, y0=0, slope=0):
        Boundary.__init__(self, type, axis, x0, y0, slope)
        return
        
    def crossed_boundary(self, pos):
        if self._slope:
            crossed = pos[self._axis] < determine_Y(pos[self._axis])
        else:
            crossed = pos[self._axis] < self._x0 
        return crossed
        
    def bounce(self, pos, vel):
        if self._slope:
            print "not supporting reflecting sloped boundaries" 
            stop
        else:
            pos[self._axis] = self._x0 + (self._x0 - pos[self._axis])
            vel[self._axis] = - vel[self._axis]
        return pos, vel
        
    def flip(self, pos, box):
        pos[self._axis] = pos[self._axis] + 2.*box[self._axis]
        return pos

class Max(Boundary):
    def __init__(self, type, axis, x0, y0=0, slope=0):
        Boundary.__init__(self, type, axis, x0, y0, slope)
        return
        
    def crossed_boundary(self, pos):
        if self._slope:
            crossed = pos[self._axis] > determine_Y(pos[self._axis])
        else:
            crossed = pos[self._axis] > self._x0 
        return crossed
        
    def bounce(self, pos, vel):
        if self._slope:
            print "not supporting reflecting sloped boundaries" 
            stop
        else:
            pos[self._axis] = self._x0 + (self._x0 - pos[self._axis])
            vel[self._axis] = - vel[self._axis]
        return pos, vel
        
    def flip(self, pos, box):
        pos[self._axis] = pos[self._axis] - 2.*box[self._axis]
        return pos
        
class PointInPlane(Boundary):
    def __init__(self, type, axis, radius, x0, y0=0, slope=0):
        Boundary.__init__(self, type, axis, x0, y0, slope)
        self._point_R2 = radius**2
        return
        
    def crossed_boundary(self, pos):
        if self._slope:
            crossed = pos[self._axis] > determine_Y(pos[self._axis])
        else:
            crossed = pos[self._axis] > self._x0 and pos[self._axis] < self._x0 + 5.*constants.nm
        return crossed
        
    def bounce(self, pos, vel):
        # does the particle pass through the hole
        if self._point_R2 < (pos[0]**2 + pos[1]**2):
            if self._slope:
                print "not supporting reflecting sloped boundaries" 
                stop
            else:
                pos[self._axis] = self._x0 + (self._x0 - pos[self._axis])
                vel[self._axis] = - vel[self._axis]
        return pos, vel
        
