import scipy
import pickle
import copy
import os
import platform
#import pylab

import Input
import particles
import normal_distribution


def loadIndexedConfs(directory="Confs", start=None, end=None):
    print directory
    indexesAndConfs = []
    confs = loadConfs(directory, start, end)
    keys = confs.keys()
    indexFirstKey = 0

    for i in xrange(0, len(confs)):
        keyFirst = keys[i]
        indexFirstKey = i
        indexTemp = range(len(confs[keyFirst]))
        if len(indexTemp):
            break
    indexesAndConfs.append(_makeIndexDict(indexTemp, confs[keyFirst]))	
    
    for i in xrange(indexFirstKey, len(confs)-1):
        indexTemp = _newIndexes(confs[keys[i+1]], confs[keys[i]], indexTemp)
        indexesAndConfs.append(_makeIndexDict(indexTemp, confs[keys[i+1]]))
    return indexesAndConfs


def _particlesAdded(confNew, confOld):
    added = 0
    mightBeAdded = (confNew[:,2] < -4.99e-7)
    
    if len(mightBeAdded):
        if mightBeAdded[-1]:
            for j in xrange(-1, -len(confNew), -1):

                diff = ((confNew[j][0]-confOld[j+added][0])**2+
                        (confNew[j][1]-confOld[j+added][1])**2+
                        (confNew[j][2]-confOld[j+added][2])**2)

                if 1e-17 < diff:
                    added = added + 1
                else:
                    break
                    
                if not mightBeAdded[j-1]:
                    break

    return added


def _particlesAbsorbed(confNew, confOld):
    absorbed = []
    mightBeAbsorbed = 4.e-7 < confOld[:,2]
    for j in xrange(len(confOld)):
        if mightBeAbsorbed[j]:
            match = 0
            j1 = j-len(absorbed)
            try :
                for i in xrange(j1, -1, -1):
                    diff = ((confOld[j][0]-confNew[i][0])**2+
                            (confOld[j][1]-confNew[i][1])**2+
                            (confOld[j][2]-confNew[i][2])**2)
                    # might be too strict!
                    if diff < 1e-16:
                        match = 1
#						print diff
                        break
            except IndexError:
                match = 0
                
            if match == 0:
                absorbed.append(j)
    return absorbed


def _newIndexes(confNew, confOld, indexOld):
    indexNew = []
    
    absorbed = _particlesAbsorbed(confNew, confOld)
    for i in xrange(len(indexOld)):
        keep = 1
        for j in absorbed:
            if i == j:
                keep = 0
                break
        if keep:
            indexNew.append(indexOld[i])
    added = _particlesAdded(confNew, confOld)

    if added:
        indexLast = indexNew[-1]+1
        for i in xrange(0, added):
            indexNew.append(i + indexLast)
    return indexNew


def _makeIndexDict(index, conf):
    dict = {}
    
    if len(index) != len(conf):
        print "The number of index and the number of atoms in conf do not match!"
        print "Nr of indexes: "+str(len(index))
        print "Nr of atoms: "+str(len(conf))
        stop
        
    for i in xrange(len(conf)):
        dict[index[i]] = conf[i]
    return dict



def followParticleR(keyPassed, zWindow, indexesAndConfs):
    posInPlane = []

    for i in xrange(len(indexesAndConfs)):
        try:
            indexesAndConfs[i][keyPassed]
            indexFirstConf = i
            break
        except KeyError:
            pass

    for i in xrange(indexFirstConf, len(indexesAndConfs)):
        print i
        try:
            posKeyPassed = indexesAndConfs[i][keyPassed]
            keysConf = indexesAndConfs[i].keys()

            for key in keysConf:
                posKey = indexesAndConfs[i][key]
                if abs(posKey[2]-posKeyPassed[2]) < zWindow:
                    posInPlane.append(posKey)
        except KeyError:
            break

    return scipy.array(posInPlane)

def followParticleV(keyPassed, zWindow, indexesAndConfs):
    velInPlane = []
    posInPlaneLastIter = {}
    
    for i in xrange(len(indexesAndConfs)):
        try:
            indexesAndConfs[i][keyPassed]
            indexFirstConf = i
            break
        except KeyError:
            pass

    for i in xrange(indexFirstConf, len(indexesAndConfs)):
        print i
        posInPlaneIter = {}
        try:
            posInPlaneIter[keyPassed] = indexesAndConfs[i][keyPassed]
            posKeyPassed = indexesAndConfs[i][keyPassed]
            keysConf = indexesAndConfs[i].keys()

            for key in keysConf:
                posKey = indexesAndConfs[i][key]
                if abs(posKey[2]-posKeyPassed[2]) < zWindow:
                    posInPlaneIter[key] = posKey

            # determine difference in positions
            keysLastIter = posInPlaneLastIter.keys()
            keysIter = posInPlaneIter.keys()

            for key1 in keysLastIter:
                for key2 in keysIter:
                    if key1 == key2:
#						vel = ((posInPlaneLastIter[key1][0]-posInPlaneIter[key1][0])**2 +
#							   (posInPlaneLastIter[key1][1]-posInPlaneIter[key1][1])**2)
                        vel = posInPlaneLastIter[key1][2]
                        velInPlane.append([posInPlaneLastIter[key1][0], posInPlaneLastIter[key1][1], vel])
            posInPlaneLastIter = posInPlaneIter

        except KeyError:
            break
    return scipy.array(velInPlane)

    
def determineAbsorbtionProfile(indexesAndConfs):
    absorbtionProfile = {}
    absorbtionProfile['space'] = []
    absorbtionProfile['time'] = [0] 
    keysLast = indexesAndConfs[0].keys()
    for i in xrange(1, len(indexesAndConfs)):
        keysNew = indexesAndConfs[i].keys()
        absorbed = []
        for keyLast in keysLast:
            match = 0
            for keyNew in keysNew:
                if keyNew == keyLast:
                    match = 1
                    break
            if match == 0:
                absorbed.append(keyLast)
        for j in absorbed:
            absorbtionProfile['space'].append(indexesAndConfs[i-1][j])
#		print absorbtionProfile['time'][-1]
        absorbtionProfile['time'].append(absorbtionProfile['time'][-1]+len(absorbed))
        keysLast = keysNew
    
    return absorbtionProfile
    
def determineEmittance(emittanceComponent, indexesAndConfs):
    if emittanceComponent == 'x' or emittanceComponent == 'X' or emittanceComponent == 0:
        emittanceComponent = 0
    if emittanceComponent == 'y' or emittanceComponent == 'Y' or emittanceComponent == 1:
        emittanceComponent = 1
    if emittanceComponent == 'z' or emittanceComponent == 'Z' or emittanceComponent == 2:
        emittanceComponent = 2
    if emittanceComponent == 'r' or emittanceComponent == 'R' or emittanceComponent == 3:
        emittanceComponent = 3

    emittancePos = 0
    emittanceVel = 0
    emittancePosVel = 0
    nrParticles = 0
        
    keysLast = indexesAndConfs[0].keys()
    for i in xrange(1, len(indexesAndConfs)):
        keysNew = indexesAndConfs[i].keys()
        absorbed = []
        for keyLast in keysLast:
            match = 0
            for keyNew in keysNew:
                if keyNew == keyLast:
                    match = 1
                    break
            if match == 0:
                absorbed.append(keyLast)

        for j in absorbed:
            pos = indexesAndConfs[i-1][j]
            vel = (indexesAndConfs[i-1][j]-indexesAndConfs[i-2][j])/(Input.timeStep*Input.sampling)
            
            if emittanceComponent == 3:
                posComponent = scipy.sqrt(pos[0]**2 + pos[1]**2) 
                velComponent = scipy.sqrt(vel[0]**2 + vel[1]**2) / vel[2]
            else:
                posComponent = pos[emittanceComponent]
                velComponent = vel[emittanceComponent]/vel[2]
            
            emittancePos = emittancePos + posComponent**2
            emittanceVel = emittanceVel + velComponent**2
            emittancePosVel = emittancePosVel + posComponent*velComponent	

            nrParticles = nrParticles+1
        keysLast = keysNew
    
    emittancePos = emittancePos/float(nrParticles)
    emittanceVel = emittanceVel/float(nrParticles)
    emittancePosVel = (emittancePosVel/float(nrParticles))**2
#	print emittancePos
#	print emittanceVel
#	print emittancePosVel
    
    return scipy.sqrt(emittancePos*emittanceVel-emittancePosVel)
    

def determineVelInPlane(planeCoor, window, axis, indexesAndConfs):
    if axis == 'x' or axis == 'X' or axis == 0:
        axisPlane = 0
        axisX = 1
        axisY = 2
    if axis == 'y' or axis == 'Y' or axis == 1:
        axisPlane = 1
        axisX = 0
        axisY = 2
    if axis == 'z' or axis == 'Z' or axis == 2:
        axisPlane = 2
        axisX = 0
        axisY = 1
        
    velInPlane = []
    posInPlaneLastIter = {}
    
    for i in xrange(len(indexesAndConfs)):
        print i
        posInPlaneIter = {}
        keysConf = indexesAndConfs[i].keys()

        for key in keysConf:
            posKey = indexesAndConfs[i][key]
            if abs(posKey[axisPlane]-planeCoor) < window:
                posInPlaneIter[key] = posKey

        # determine difference in positions
        keysLastIter = posInPlaneLastIter.keys()
        keysIter = posInPlaneIter.keys()

        for key1 in keysLastIter:
            for key2 in keysIter:
                if key1 == key2:
                    vel = ((posInPlaneLastIter[key1][axisX]-posInPlaneIter[key1][axisX])**2 +
                           (posInPlaneLastIter[key1][axisY]-posInPlaneIter[key1][axisY])**2)# + 
#						   (posInPlaneLastIter[key1][axisPlane]-posInPlaneIter[key1][axisPlane])**2)
#					vel = (posInPlaneLastIter[key1][2]-posInPlaneIter[key1][2])**2
#					vel = (posInPlaneLastIter[key1][axisPlane]-posInPlaneIter[key1][axisPlane])**2
                    velInPlane.append([posInPlaneLastIter[key1][axisX], posInPlaneLastIter[key1][axisY], vel])

        posInPlaneLastIter = posInPlaneIter

#		except KeyError:
#			break

    return scipy.array(velInPlane)
    

def determineDensityAboveEmitter(window, indexesAndConfs):
    densityAll = []	
    for i in xrange(len(indexesAndConfs)):
        density = 0
        rEmitter2 = Input.rDeposit**2
        particleInWindow = {}
        keysConf = indexesAndConfs[i].keys()

        for key in keysConf:
            posKey = indexesAndConfs[i][key]
            if abs(posKey[2]-(-Input.zDeposit)) < window:
                r2 = posKey[0]**2+posKey[1]**2
                if r2 < rEmitter2:
                    density = density + 1
        densityAll.append(density/(scipy.pi*rEmitter2*window)*Input.coulombPerElectron)
    return densityAll


def getConf(index, indexesAndConfs):
    conf = indexesAndConfs[index]
    particles = particles.Particles(len(conf))
    pos = []
    for key in conf.keys():
        pos.append(conf[key])
    particles.setR(scipy.array(pos))
    return particles


def determineEFielndInPoint(pos, particles):
    particlesTest = copy.copy(particles)
    particlesTest.add(pos)
    indexLast = len(particlesTest)-1
    return particlesTest.getSingleF(indexLast)/1.6021e-19

    
def determineEFieldInPlane(resolution, particles, axis=0):
    meshEField = scipy.zeros((resolution, resolution),"f")
    dx = 2.*Input.x/resolution/1.
    dy = 2.*Input.y/resolution/1.
    
    for i in xrange(0,resolution):
        for j in xrange(0,resolution):
            pos = [-Input.x/1.+i*dx, -Input.y/1.+j*dy, 0.]#Input.z+1e-8]
            field = determineEFielndInPoint(pos, particles)[2]
            field = scipy.log10(field)
            if 6.5 < field:
                field = 6.5
            if field < 1.5:
                field = 1.5
                
            meshEField[i,j] = field
    meshEField[0,0] = 7
    meshEField[resolution-1,resolution-1] = 1.5

    return meshEField


def determineRawAbsorptionTimeProfile(allParticlesAbsorbed, resolution=0):
    lastInterval = 0
    timeProfile = []

    if resolution:
        interval = len(allParticlesAbsorbed)/resolution	
        for i in xrange(0, resolution):
            timeProfile.append(allParticlesAbsorbed[int(i*interval)]-lastInterval)
            lastInterval = allParticlesAbsorbed[int(i*interval)]
    else:
        for i in xrange(0, len(allParticlesAbsorbed)):
            timeProfile.append(allParticlesAbsorbed[i]-lastInterval)
            lastInterval = allParticlesAbsorbed[i]
            
    return scipy.array(timeProfile)


def determineSmoothedAbsorptionTimeProfile(indexesAndConfs, interval):
    timeProfile = []
    absorbed = determineAbsorbtionProfile(indexesAndConfs)
#	absorbed = data[3]
#	allParticlesAbsorbed = determineRawAbsorptionTimeProfile(absorbed['time'], 1000)
    allParticlesAbsorbed = determineRawAbsorptionTimeProfile(absorbed['time'])
    
    for i in xrange(interval, len(allParticlesAbsorbed)-interval):
        averageInterval = 0
        for ii in xrange(-interval, interval+1):
            if allParticlesAbsorbed[i+ii]:
                averageInterval = averageInterval+allParticlesAbsorbed[i+ii]*normal_distribution.normal(ii,0.,interval/2.)
        timeProfile.append(averageInterval)
                
    return scipy.array(timeProfile)		

    
def determinePulsSizeAndEmittance(dirListConf):
    pulsSizeAndEmittance = []
    for dirConf in dirListConf:
        c = loadIndexedConfs(dirConf)
        ap = determineAbsorbtionProfile(c)
        e = determineEmittance(3,c)
        pulsSizeAndEmittance.append([ap['time'][-1], e])
    return pulsSizeAndEmittance
    
def savePulsSizeAndEmittance(dirListConf):
    f0 = open("puls.txt","w")
    f1 = open("emittance.txt","w")
    for dir in dirListConf:
        data = determinePulsSizeAndEmittance([dir])
#		print data[0][0]
#		print data
        f0.write(str(data[0][0]) + "\n")
        f1.write(str(data[0][1]) + "\n")
    f0.close()
    f1.close()
    
    f2 = open("intensity.txt","w")
    for d in dirListConf:
        intensity = 0
        for i in range(0,len(d)):
            if d[i] == "_":
                intensity = d[i+1:]
        f2.write(str(intensity) + "\n")
    f2.close()
    return
    

def saveSmoothedAbsorptionTimeProfile(dirListConf):
    for dir in dirListConf:
        c = loadIndexedConfs(dir)
        timeProfile = determineSmoothedAbsorptionTimeProfile(c,200)

        f = open(dir+"timeProfile.txt","w")
        for t in timeProfile:
            f.write(str(t) + "\n")
        f.close()
    return


def determineHistogramZ(conf, bins):
    z = []
    for key in conf.keys():
        z.append(conf[key][2])
                
    pylab.hist(z, bins, [-Input.z, Input.z], align='mid')
    pylab.xlim([-Input.z, Input.z])
    pylab.ylim([0, 50])
    return


def saveSeriesOfHistograms(confs, bins):
    i = 0
    for j in xrange(0, len(confs), 10):
        pylab.clf()
        determineHistogramZ(confs[j], bins)
        if i < 10:
            pylab.savefig("hist/Fig000%i.png"%i)
        elif i< 100:
            pylab.savefig("hist/Fig00%i.png"%i)
        elif i< 1000:
            pylab.savefig("hist/Fig0%i.png"%i)
        else:
            pylab.savefig("hist/Fig%i.png"%i)
        i = i + 1
    return	

def determineSmoothedZ(conf, ii):
    zAll = []
    for key in conf.keys():
        zAll.append(conf[key][2])

    zSmoothAll = []
    dz = 2.*Input.z/float(ii)
    
    for i in xrange(ii):
        zSmooth = 0
        zCurrentIter = -Input.z + i * dz
        for z in zAll:
            zSmooth = zSmooth + normal_distribution.normal(z - zCurrentIter, 0, 10.*Input.nm)
        zSmoothAll.append(zSmooth)
    return scipy.array(zSmoothAll)		

    
def saveSeriesOfSmoothed(confs, ii):
    i = 0
    for j in xrange(0, len(confs), 10):
        pylab.clf()
        dataSmooth = determineSmoothedZ(confs[j], ii)
        pylab.plot(dataSmooth)
        pylab.ylim([0,5e8])
        pylab.text(450,4.5e8,str(len(confs[j].keys())))
        if len(confs[j].keys()) == 0:
            break
        if i < 10:
            pylab.savefig("smooth/Fig000%i.png"%i)
        elif i< 100:
            pylab.savefig("smooth/Fig00%i.png"%i)
        elif i< 1000:
            pylab.savefig("smooth/Fig0%i.png"%i)
        else:
            pylab.savefig("smooth/Fig%i.png"%i)
        i = i + 1
    return
        
def loadConfigurationDirectories(path="."):
    dirListAll = os.listdir(path)
    dirListConfTemp = []
    for fname in dirListAll:
        if (os.path.isdir(fname) and fname[0:15] == "Configurations_" and fname[15:].isalnum()):
            dirListConfTemp.append(fname[15:])
    dirListConfTemp = scipy.array(dirListConfTemp,"i")
    dirListConfTemp.sort()
    dirListConf = []
    for dirConfTemp in dirListConfTemp:
        dirListConf.append("Configurations_"+str(dirConfTemp))
    return dirListConf
    
    
def loadDataAll():
    pulses = []
    emittances = []
    
    f = open("dataAll1.txt","r")
    lines = f.readlines()
    f.close()
    
    for line in lines:
        for i in xrange(len(line)):
            if "[" == line[i] and line[i+1].isalnum():
                startPulse = i+1
            if "," == line[i] and line[i+1] == " " and line[i+2].isalnum():
                endPulse = i
                startEmittance = i+2
            if line[i].isalnum() and "]" == line[i+1]:
                endEmittance = i+1
        pulses.append(int(line[startPulse:endPulse]))
        emittances.append(float(line[startEmittance:endEmittance]))
    return pulses, emittances
    
    
def loadConfs(directory="Confs", start=None, end=None):
    if platform.system() == 'Windows':
        directory = directory + "\\"
    elif platform.system() == 'Linux':
        directory = directory + "/"

    confs = {}
    if start == None:
        start = 0
    if end == None:
        end = 1000000
    try:
        for i in range(start,end):
            file = open(directory+str(i)+".pickle","r")
            confs[i+start] = pickle.load(file)
    except IOError:
        pass
    return confs
