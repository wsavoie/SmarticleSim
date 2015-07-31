import numpy
import scipy
import sys, os, errno, platform,shutil
import winsound,time, datetime
from time import gmtime, strftime
from datetime import date
import ctypes

compName = platform.node(); 

def chooseDir(name): #add your computer name in the dictionary here and the corresponding location where the program output will be placed
     return{
     'PHYS32240':'D:\\SimResults\\Chrono\\SmarticleU\\Results\\',
     'WS':'C:\\SimResults\\Chrono\\SmarticleU\\Results\\'}[name]
def chooseRunLoc(name): #add your computer name in the dictionary here and the corresponding location of the exe file
     return{
     'PHYS32240':'"D:\ChronoCode\chronoPkgs\SmarticlesBuild\Release\SmarticlesSystem.exe"',
     'WS':'"C:\ChronoCode\chronoPkgs\SmarticlesBuild\Release\SmarticlesSystem.exe"'}[name]
     
def chooseParLoc(name): #add your computer name in the dictionary here and the corresponding location of the exe file
    return{
    'PHYS32240':'D:\ChronoCode\chronoPkgs\Smarticles\pythonrunscripts\simRunPars.txt',
    'WS':'"C:\ChronoCode\chronoPkgs\SmarticlesBuild\Release\SmarticlesSystem.exe"'}[name]
def makePath(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise Exception('path error problem')
            
def getFileNum():
    x =os.path.basename(__file__)
    return{'smarticle_run.py':0,'smarticle_run2.py':1,'smarticle_run3.py':2
    }[x]
def getPars():
    with open(chooseParLoc(compName), 'r+') as f:
        read_data = f.readlines()
    f.closed
    fnum = getFileNum() #returns number based on current file
    
    dt=read_data[0]
    dt=float(dt)
    
    angle= read_data[fnum+2]
    angle=[int(x) for x in angle.split('\t')]
    
    lw = read_data[fnum+6]
    lw = [float(x) for x in lw.split('\t')]
    
    numlayer = read_data[fnum+10]
    numlayer = [int(x) for x in numlayer.split('\t')]
    
    gamma = read_data[fnum+14]
    gamma = [float(x) for x in gamma.split('\t')]
    
    read = read_data[fnum+18]
    read = [int(x) for x in read.split('\t')]
    
    return [dt, angle, lw, numlayer, gamma, read]
def runSim():
    cores = 1
    sliding_its = 55
    bilateral_its= 55
    time.mktime
    [dT, angle, lw, numlayers, gamma, read] = getPars()
    
    #set depVar to something just in case first 2 
    if (len(angle)>len(lw)):
        depVar=angle
    elif (len(lw)>len(angle)):
        depVar=lw
    elif (read[0]==1): #this param takes priority
        depVar=gamma
    else:#this may change if I add another type of parameter
        depVar=gamma
        
        
    if (depVar==angle):
        lw = lw*len(depVar)
        numlayers = numlayers*len(depVar)
        gamma = gamma*len(depVar)
    elif (depVar==lw): #will define numlayers explicitly in this case
        angle = angle*len(depVar)
        gamma = gamma*len(depVar)
    elif (depVar==gamma):
        angle = angle*len(depVar)
        lw = lw*len(depVar)
        numlayers = numlayers*len(depVar)
    #make file
    print read
    for i in range(0,len(depVar)):
        d=date.fromtimestamp(time.time())
        t = d.strftime("%Y%m%d")
        dirn = "%s lw=%g ang=%g gamma=%g "%(t, lw[i], angle[i],gamma[i])+read[0]*"r"
        dirpath = chooseDir(compName)+dirn
        print dirpath
        makePath(dirpath)
        shutil.copyfile(chooseParLoc(compName)+"\..\smarticles.csv", dirpath+"\smarticles.csv")
        os.chdir(dirpath)
        tBegin = time.time()
        
        simT = time.time();          
        x= "%s %f %g %g %g %g %g %g %g %g"%(fileloc,lw[i], dT,numlayers[i],angle[i], gamma[i], read[0], cores,sliding_its,bilateral_its)
        title= "%g %g %g %g %g %g %g %g %g %g"%(getFileNum()+1,lw[i], dT,numlayers[i],angle[i],gamma[i], read[0],cores,sliding_its,bilateral_its)
        ctypes.windll.kernel32.SetConsoleTitleA(title)
        print 'hi'
        print x
        print 'hi'
        os.system(x)

        endSimt = time.time()-simT
        print "######################"
        print strftime("%H:%M:%S",gmtime())
        print "######################"
        winsound.Beep(1000,300)
        tElapsed = time.time()-tBegin
        os.chdir("..\\..\\..")
        
        # renaming folder upon completion doesnt work
        ofpath = dirpath
        fpath = dirpath
        count = 1
        while True:
            try:
                # split=os.path.splitext(fpath)
                # fpath= split[0]+'b'+split[1]
                fpath = ofpath+' '+str(count)
                os.rename(ofpath,fpath)
                break
            except (NameError,WindowsError):
                print 'trying rename again'
                count+=1
fileloc=chooseRunLoc(compName)
runSim()

#if__name__=="__main__":
#	sys.exit(main())


