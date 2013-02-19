import scipy
import pylab

import setup
reload(setup)

import analyze_tools
reload(analyze_tools)

import constants
reload(constants)

plotFrame = 5*constants.nm

def save(posAll, indexStart=0, indexEnd=None):
    if indexEnd is None:
        indexEnd = len(posAll)
    for i in xrange(indexStart,indexEnd):
        pylab.clf()
        plotSingleXZ(posAll,i )
        
        print str(i) + " of "+str(indexEnd) 

        if i < 10:
            pylab.savefig("fig/Fig000%i.png"%i)
        elif i< 100:
            pylab.savefig("fig/Fig00%i.png"%i)
        elif i< 1000:
            pylab.savefig("fig/Fig0%i.png"%i)
        else:
            pylab.savefig("fig/Fig%i.png"%i)
    return


def plotSingleXY(posAll, i):
    ax = pylab.fill([-setup.x-plotFrame,setup.x+plotFrame,setup.x+plotFrame,-setup.x-plotFrame],[-setup.y-plotFrame,-setup.y-plotFrame,setup.y+plotFrame,setup.y+plotFrame],'r')
    bx = pylab.fill([-setup.x,setup.x,setup.x,-setup.x],[-setup.y,-setup.y,setup.y,setup.y],'w')
    pylab.plot([posAll[i][:,0], posAll[i][:,0]], 
               [posAll[i][:,1], posAll[i][:,1]],
               'k.', markersize=5.)
    return

def plotSingleXZ(posAll, i):
    ax = pylab.fill([-setup.x-plotFrame,setup.x+plotFrame,setup.x+plotFrame,-setup.x-plotFrame],[-setup.z-plotFrame,-setup.z-plotFrame,setup.z+plotFrame,setup.z+plotFrame],'r')
    bx = pylab.fill([-setup.x,setup.x,setup.x,-setup.x],[-setup.z,-setup.z,setup.z,setup.z],'w')
    pylab.plot([posAll[i][:,0], posAll[i][:,0]], 
               [posAll[i][:,2], posAll[i][:,2]],
               'k.', markersize=5.)
    return

def plotSingleYZ(posAll, i):
    ax = pylab.fill([-setup.y-plotFrame,setup.y+plotFrame,setup.y+plotFrame,-setup.y-plotFrame],[-setup.z-plotFrame,-setup.z-plotFrame,setup.z+plotFrame,setup.z+plotFrame],'r')
    bx = pylab.fill([-setup.y,setup.y,setup.y,-setup.y],[-setup.z,-setup.z,setup.z,setup.z],'w')
    pylab.plot([posAll[i][:,1], posAll[i][:,1]], 
               [posAll[i][:,2], posAll[i][:,2]],
               'r.', markersize=5.)
    return


def plotAbsorbed(totalParticlesAbsorbed):
    totalParticlesAbsorbed = scipy.array(totalParticlesAbsorbed)
    ax = pylab.fill([-setup.x-plotFrame,setup.x+plotFrame,setup.x+plotFrame,-setup.x-plotFrame],[-setup.y-plotFrame,-setup.y-plotFrame,setup.y+plotFrame,setup.y+plotFrame],'r')
    bx = pylab.fill([-setup.x,setup.x,setup.x,-setup.x],[-setup.y,-setup.y,setup.y,setup.y],'w')
    pylab.plot([totalParticlesAbsorbed[:,0], totalParticlesAbsorbed[:,0]], 
               [totalParticlesAbsorbed[:,1], totalParticlesAbsorbed[:,1]],
               'k.', markersize=5.)
    return
	
def plotSampling(data):
    x = scipy.array(range(0, len(data)))
#	print x
    x = x*setup.timeStep*setup.sampling
    pylab.plot(x, data)
    pylab.xlabel('Time')
    return

def flipThrough(posAll, startIndex=0, stepSize=1):
    i = 0
    pylab.clf()
    while 1:
        plotSingleXZ(posAll, i+startIndex)
        i = i+stepSize
        setup = raw_setup()
        pylab.clf()
    return
    
def plotPoints(data):
    a = scipy.array(data)
    s = scipy.zeros(len(a))+0.01
    pylab.scatter(a[:,0], a[:,1], s)
    pylab.axis('equal')
    return

def plotTimeEvolutionEField(indexesAndConfs, resolution=50):
    for i in xrange(0,len(indexesAndConfs)):
        print i
        p = analyze_tools.recreateConf(i,indexesAndConfs)
        m = analyze_tools.determineEFieldInPlane(resolution,p)
        pylab.clf()
        pylab.pcolormesh(m)
        pylab.colorbar()
        if i < 10:
            pylab.savefig("00"+str(i)+"E_field.png")
        elif i < 100:
            pylab.savefig("0"+str(i)+"E_field.png")
        else:
            pylab.savefig(str(i)+"E_field.png")
