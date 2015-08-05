import sys, os, errno, platform,shutil
import time, datetime
from time import gmtime, strftime
from datetime import date
import ctypes

compName = platform.node()

def chooseDir(name): #add your computer name in the dictionary here and the corresponding location where the program output will be placed
     return{
     'PHYS32240':'D:\\SimResults\\Chrono\\SmarticleU\\Results\\',       #lab computer
     'PHYS32164':'D:\\LabDropbox\\Dropbox\\WillSmarticles\\Results\\',    #goldman2
     'WS':'C:\\SimResults\\Chrono\\SmarticleU\\Results\\',              #laptop
     'euler.wacc.wisc.edu':'/home/wsavoie/SmarticleResults/'            #Wisc server
     }[name]
def chooseRunLoc(name): #add your computer name in the dictionary here and the corresponding location of the exe file
     return{
     'PHYS32240':'"D:\ChronoCode\chronoPkgs\SmarticlesBuild\Release\SmarticlesSystem.exe"',
     'PHYS32164':'D:\\LabDropbox\\Dropbox\\WillSmarticles\\Release\\SmarticlesSystem.exe',
     'WS':'"C:\ChronoCode\chronoPkgs\SmarticlesBuild\Release\SmarticlesSystem.exe"',
     'euler.wacc.wisc.edu':'/home/wsavoie/ChronoSrc/SmarticlesBuild/SmarticlesSystem'}[name]
     
def chooseParLoc(name): #add your computer name in the dictionary here and the corresponding location of the exe file
    return{
    'PHYS32240':'D:\ChronoCode\chronoPkgs\Smarticles\pythonrunscripts\simRunPars.txt',
    'PHYS32164':'D:\LabDropbox\Dropbox\WillSmarticles\Smarticles\pythonrunscripts\simRunPars.txt',
    'WS':'"C:\ChronoCode\chronoPkgs\SmarticlesBuild\Release\SmarticlesSystem.exe"',
    'euler.wacc.wisc.edu':'/home/wsavoie/ChronoSrc/Smarticles/pythonrunscripts/simRunPars.txt'}[name]
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
    
    angle1= read_data[fnum+2]
    angle1=[int(x) for x in angle1.split('\t')]
    
    angle2= read_data[fnum+6]
    angle2=[int(x) for x in angle2.split('\t')]
    
    lw = read_data[fnum+10]
    lw = [float(x) for x in lw.split('\t')]
    
    numlayer = read_data[fnum+14]
    numlayer = [int(x) for x in numlayer.split('\t')]
    
    gamma = read_data[fnum+18]
    gamma = [float(x) for x in gamma.split('\t')]
    
    read = read_data[fnum+22]
    read = [int(x) for x in read.split('\t')]
    
    return [dt, angle1,angle2, lw, numlayer, gamma, read]
def runSim():
    cores = 1
    sliding_its = 55
    bilateral_its= 55
    time.mktime
    [dT, angle1,angle2, lw, numlayers, gamma, read] = getPars()
    
    #break if len(angle1) != len(angle2)
    if(len(angle1)!=len(angle2)):
        sys.exit("error in simRunPars, angle inputs are not equal in length")
    #set depVar to something just in case first 2 
    if (len(angle1)>len(lw)):
        depVar=angle1
    elif (len(lw)>len(angle1)):
        depVar=lw
    elif (read[0]==1): #this param takes priority
        depVar=gamma
    else:#this may change if I add another type of parameter
        depVar=gamma
        
        
    if (depVar==angle1):
        lw = lw*len(depVar)
        numlayers = numlayers*len(depVar)
        gamma = gamma*len(depVar)
    elif (depVar==lw): #will define numlayers explicitly in this case
        angle1 = angle1*len(depVar)
        angle2 = angle2*len(depVar)
        gamma = gamma*len(depVar)
    elif (depVar==gamma):
        angle1 = angle1*len(depVar)
        angle2 = angle2*len(depVar)
        lw = lw*len(depVar)
        numlayers = numlayers*len(depVar)
    #make file
    print read
    for i in range(0,len(depVar)):
        d=date.fromtimestamp(time.time())
        t = d.strftime("%Y%m%d")
        dirn = "%s lw=%g ang1=%g ang2=%g g=%g "%(t, lw[i], angle1[i],angle2[i],gamma[i])+read[0]*"r"
        dirpath = chooseDir(compName)+dirn
        print dirpath
        makePath(dirpath)
        if read[0]==1:
            shutil.copyfile(chooseParLoc(compName)+"\..\smarticles.csv", dirpath+"\smarticles.csv")
        os.chdir(dirpath)
        tBegin = time.time()
        
        simT = time.time();          
        x= "%s %f %g %g %g %g %g %g %g %g %g"%(fileloc,lw[i], dT,numlayers[i],angle1[i],angle2[i], gamma[i], read[0], cores,sliding_its,bilateral_its)
        x2= "%f %g %g %g %g %g %g %g %g %g"%(lw[i], dT,numlayers[i],angle1[i],angle2[i], gamma[i], read[0], cores,sliding_its,bilateral_its)
        title= "%g %g %g %g %g %g %g %g %g %g %g"%(getFileNum()+1,lw[i], dT,numlayers[i],angle1[i],angle2[i],gamma[i], read[0],cores,sliding_its,bilateral_its)
        # ctypes.windll.kernel32.SetConsoleTitleA(title)
        print 'hi'
        print x
        print 'hi'
        if sys.platform != 'win32':
            os.system('qsub /home/wsavoie/ChronoSrc/Smarticles/pythonrunscripts/bash_Smarticle.sh -F "'+x2+'"')
        else:
            os.system(x)
        
        endSimt = time.time()-simT
        print "######################"
        print strftime("%H:%M:%S",gmtime())
        print "######################"
        # winsound.Beep(1000,300)
        tElapsed = time.time()-tBegin
        os.chdir("../../")
        
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


