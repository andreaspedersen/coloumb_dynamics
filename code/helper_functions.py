import numpy

def normal_distribution(x, mhy, sigma):
    return 1./(sigma*numpy.sqrt(2.*numpy.pi))*numpy.exp(-((x-mhy)**2)/(2.*sigma**2))
    
def inverse_normal_distribution(x, mhy, sigma):
    return (sigma*numpy.sqrt(2.*numpy.pi))/numpy.exp(-((x-mhy)**2)/(2.*sigma**2))

def make_gaussian_distribution(tInterval, nrParticles):
    xMin = -2.5
    xMax =  2.5
    mhy = 0.0
    sigma = 1.0
    dx = (xMax - xMin) / float(nrParticles)
    emissionRate = numpy.zeros(nrParticles,'f')
    for i in xrange(0, nrParticles):
        emissionRate[i] = inverse_normal_distribution(xMin+i*dx, mhy, sigma)
    emissionRate = tInterval*emissionRate/sum(emissionRate)
    return emissionRate
    
def make_square_distribution(tInterval, nrParticles):
    dt = (tInterval) / float(nrParticles)
    emissionRate = numpy.zeros(nrParticles,'f')
    for i in xrange(0, nrParticles):
        emissionRate[i] = dt
    return emissionRate
