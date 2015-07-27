import numpy
import scipy
import sys, os, errno, platform
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
    
    return [dt, angle, lw, numlayer]
def runSim():
    cores = 1
    sliding_its = 55
    bilateral_its= 55
    time.mktime
    [dT, angle, lw, numlayers] = getPars()
    
    if (len(angle)>len(lw)):
        depVar=angle
    else:
        depVar=lw
        
    if (depVar==angle):
        lw = lw*len(angle)
        numlayers = numlayers*len(angle)
    elif (depVar==lw):
        angle = angle*len(depVar)
    #make file
    
    for i in range(0,len(depVar)):
        d=date.fromtimestamp(time.time())
        t = d.strftime("%Y%m%d")
        dirn = "%s lw=%g ang=%g"%(t, lw[i], angle[i])
        dirpath = chooseDir(compName)+dirn
        print dirpath
        makePath(dirpath)
        os.chdir(dirpath)
        tBegin = time.time()
        
        simT = time.time();          
        x= "%s %.2f %g %g %g %g %g %g"%(fileloc,lw[i], dT,numlayers[i],angle[i],cores,sliding_its,bilateral_its)
        title= "%g %g %g %g %g %g %g %g"%(getFileNum()+1,lw[i], dT,numlayers[i],angle[i],cores,sliding_its,bilateral_its)
        ctypes.windll.kernel32.SetConsoleTitleA(title)
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


