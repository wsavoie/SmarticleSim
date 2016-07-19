set LW=.8
set DT=.00025
set NUMLAYER=2
REM 2=R&W 1=R 0='' -1=W
set READFILE=0
set ACTIVEPCT=1
set ANGLE1=45
set ANGLE2=45
set BOXANG=-30
set NUMPERLAYER=20
set STRESSACTIVATE=0.3
set SCREENSHOT=0
start "SmarticlesSystem_Dynamic.exe" /high "D:\ChronoCode\chronoPkgs\SmarticlesBuild\Release\SmarticlesSystem_Dynamic.exe" %LW% %DT% %NUMLAYER% %READFILE% %ACTIVEPCT% %ANGLE1% %ANGLE2% %BOXANG% %NUMPERLAYER% %STRESSACTIVATE% %SCREENSHOT%