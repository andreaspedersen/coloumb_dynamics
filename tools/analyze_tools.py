import scipy
import pickle
import copy
import platform

import setup
import particles
import constants

def saveList(dataList, filename="dataList.txt"):
    f = open(filename, "w")
    for dataElement in dataList:
        f.write(str(dataElement) + "\n")
    f.close()
    return
    
def normal(x, mhy, sigma):
    return 1./(sigma*scipy.sqrt(2.*scipy.pi))*scipy.exp(-((x-mhy)**2)/(2.*sigma**2))

#def determineCurrentDensity(totalParticlesAbsorbed):
#    current = []
#    steps = xrange(1,len(totalParticlesAbsorbed)+1)
#    current = scipy.zeros(len(totalParticlesAbsorbed))
#    for i in xrange(0,len(steps)):
#        current[i] = (totalParticlesAbsorbed[i])/(2.*setup.x*2.*setup.y)/(setup.timeStep*(i+1))*constants.coulombPerElectron
#    print current[-1]
#    print 2.3e-6*(setup.V_z_constant)**(3./2.)/(2.*setup.z)**2
#    print current[-1]/(2.3e-6*(setup.V_z_constant)**(3./2.)/(2.*setup.z)**2)
#    return current

def determineAbsorptionProfilePlane(allParticlesAbsorbed, resolution):
    meshAbsorbed = scipy.zeros((resolution, resolution),"f")
    dx = 2.*setup.x/resolution
    dy = 2.*setup.y/resolution
    
    for particleAbsorbed in allParticlesAbsorbed:
        ix = (setup.x + particleAbsorbed[0])/dx
        iy = (setup.y + particleAbsorbed[1])/dy
        if resolution <= ix:
            ix = resolution - 1
        if resolution <= iy:
            iy = resolution - 1
        meshAbsorbed[iy,ix] = meshAbsorbed[iy,ix] + 1
        
    return meshAbsorbed

def determineAbsorptionProfileLine(allParticlesAbsorbed, resolution):
    lineAbsorbed = scipy.zeros(resolution,"f")
    dr = scipy.sqrt(setup.x**2+setup.y**2) / resolution
    
    for particleAbsorbed in allParticlesAbsorbed:
        ir = scipy.sqrt(particleAbsorbed[0]**2+particleAbsorbed[1]**2) / dr
#		# to catch particles with a radius larger than the box length
#		if resolution <= ir:
#			ir = resolution - 1
        lineAbsorbed[ir] = lineAbsorbed[ir] + 1

    for i in xrange(len(lineAbsorbed)):
        lineAbsorbed[i] = lineAbsorbed[i] / (2.*scipy.pi*(i+1.)*dr)
    return lineAbsorbed	

def smoothAbsorptionTimeProfile(data, interval):
    timeProfile = []
    absorbed = data[3]
    allParticlesAbsorbed = determineAbsorptionTimeProfile(absorbed['nr'], 1000)
    
    for i in xrange(interval, len(allParticlesAbsorbed)-interval):
        averageInterval = 0
        for ii in xrange(-interval, interval+1):
            if allParticlesAbsorbed[i+ii]:
                averageInterval = averageInterval+allParticlesAbsorbed[i+ii]*normal(ii,0.,interval/2.)
        timeProfile.append(averageInterval)
                
    return scipy.array(timeProfile)	
    
def determineAbsorptionTimeProfile(allParticlesAbsorbed, resolution):
    interval = len(allParticlesAbsorbed)/resolution
    lastInterval = 0
    timeProfile = []
    
    for i in xrange(0, resolution):
        timeProfile.append(allParticlesAbsorbed[int(i*interval)]-lastInterval)
        lastInterval = allParticlesAbsorbed[int(i*interval)]
                
    return scipy.array(timeProfile)
    
def determineAverageAbsorptionTimeProfile(dataAll, resolution):
    absorbed = dataAll[0][3]
    average = determineAbsorptionTimeProfile(absorbed['nr'], resolution)
    for i in xrange(1,len(dataAll)):
        absorbed = dataAll[i][3]
        average = average + determineAbsorptionTimeProfile(absorbed['nr'], resolution)
    average = average/float(len(dataAll))
    return average
    
def particlesAdded(confNew, confOld):
    added = 0
    mightBeAdded = (confNew[:,2] < -4.99e-7)
    
    if len(mightBeAdded):
        if mightBeAdded[-1]:
            for j in xrange(-1, -len(confNew), -1):

                diff = ((confNew[j][0]-confOld[j+added][0])**2+
                        (confNew[j][1]-confOld[j+added][1])**2+
                        (confNew[j][2]-confOld[j+added][2])**2)

                if 1e-16 < diff:
                    added = added + 1
                else:
                    break
                    
                if not mightBeAdded[j-1]:
                    break
    #	print "Added "+str(added)

    return added

def particlesAbsorbed(confNew, confOld):
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
                    if diff < 1e-13:
                        match = 1
                        break
            except IndexError:
                match = 0
                
            if match == 0:
                absorbed.append(j)
#	print "Absorbed "+str(absorbed)
    return absorbed

def determineNewIndexes(confNew, confOld, indexOld):
    indexNew = []
        
    absorbed = particlesAbsorbed(confNew, confOld)
#	print absorbed
    for i in xrange(len(indexOld)):
        keep = 1
        for j in absorbed:
            if i == j:
#				print indexOld[i]
                keep = 0
                break
        if keep:
            indexNew.append(indexOld[i])

    print len(indexNew)
    
    indexLast = indexNew[-1]+1
    added = particlesAdded(confNew, confOld)

    print added

    for i in xrange(0, added):
        indexNew.append(i + indexLast)

    print len(indexNew)
    print len(confNew)
    print len(confOld)
    print "-------"
    return indexNew
    
def makeIndexDict(index, conf):
    dict = {}
    
    for i in xrange(len(conf)):
        dict[index[i]] = conf[i]
    return dict

def savePosAll(posAll, directory="Confs", start=None, end=None):
    if platform.system() == 'Windows':
        directory = directory + "\\"
    elif platform.system() == 'Linux':
        directory = directory + "/"

    if start == None and end == None:
        interval = len(posAll)
        start = 0
    else:
        interval = end-start

    for i in xrange(interval):
        file = open(directory+str(i+start)+".pickle","w")
        pickle.dump(posAll[i+start], file)
    return

    
def loadSequence(start, end, directory="Confs"):
    if platform.system() == 'Windows':
        directory = directory + "\\"
    elif platform.system() == 'Linux':
        directory = directory + "/"

    confs = {}
    interval = end-start
    for i in xrange(interval):
        file = open(directory+str(i+start)+".pickle","r")
        confs[i+start] = pickle.load(file)
    return confs
    
def loadIndexedSequence(start, end, directory="Confs"):
    if platform.system() == 'Windows':
        directory = directory + "\\"
    elif platform.system() == 'Linux':
        directory = directory + "/"

    indexesAndConfs = []
    confs = loadSequence(start, end, directory)
    keys = confs.keys()
    indexFirstKey = 0

    for i in xrange(0, len(confs)):
        keyFirst = keys[i]
        indexFirstKey = i
        indexTemp = range(len(confs[keyFirst]))
        if len(indexTemp):
            break
    indexesAndConfs.append(makeIndexDict(indexTemp, confs[keyFirst]))	
    
    for i in xrange(indexFirstKey, len(confs)-1):
        print "------"
        print "index "+str(i)

        indexTemp = determineNewIndexes(confs[keys[i+1]], confs[keys[i]], indexTemp)

        print len(confs[keys[i+1]])
        print len(indexTemp)

        indexesAndConfs.append(makeIndexDict(indexTemp, confs[keys[i+1]]))

        print "------"

    return indexesAndConfs

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

def planeV(planeCoor, window, axis, indexesAndConfs):
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
    
def densityAboveEmitter(window, indexesAndConfs):
    densityAll = []	
    for i in xrange(len(indexesAndConfs)):
        density = 0
        rEmitter2 = setup.rDeposit**2
        particleInWindow = {}
        keysConf = indexesAndConfs[i].keys()

        for key in keysConf:
            posKey = indexesAndConfs[i][key]
            if abs(posKey[2]-(-setup.zDeposit)) < window:
                r2 = posKey[0]**2+posKey[1]**2
                if r2 < rEmitter2:
                    density = density + 1
        densityAll.append(density/(scipy.pi*rEmitter2*window)*constants.coulombPerElectron)
    return densityAll

def recreateConf(index, indexesAndConfs):
    conf = indexesAndConfs[index]
    particles = Particles.Particles(len(conf))
#	keys = conf.keys()
    pos = []
    for key in conf.keys():
        pos.append(conf[key])
#	print len(particles)
#	print len(pos)
    particles.set_R(scipy.array(pos))
    return particles

def getEFielndInPoint(pos, particles):
    particlesTest = copy.copy(particles)
    particlesTest.add(pos)
    indexLast = len(particlesTest)-1
    return particlesTest.getSingleF(indexLast)/1.6021e-19
    
def determineEFieldInPlane(resolution, particles, axis=0):
    meshEField = scipy.zeros((resolution, resolution),"f")
    dx = 2.*setup.x/resolution/1.
    dy = 2.*setup.y/resolution/1.
    
    for i in xrange(0,resolution):
        for j in xrange(0,resolution):
            pos = [-setup.x/1.+i*dx, -setup.y/1.+j*dy, 0.]#setup.z+1e-8]
            field = getEFielndInPoint(pos, particles)[2]
            field = scipy.log10(field)
            if 6.5 < field:
                field = 6.5
            if field < 1.5:
                field = 1.5
                
            meshEField[i,j] = field
    meshEField[0,0] = 7
    meshEField[resolution-1,resolution-1] = 1.5

    return meshEField
