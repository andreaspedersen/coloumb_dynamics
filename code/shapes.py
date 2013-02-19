class Line:
    def __init__(self, x0, y0, slope):
        self._x0 = x0
        self._y0 = y0
        self._slope = slope
        return

    def determine_Y(self, x):
        return self._y0 + (x-self._x0)*self._slope