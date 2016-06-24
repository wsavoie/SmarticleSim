set LW=1
set DT=.00025
set NUMLAYER=1
REM 2=R&W 1=R 0='' -1=W
set READFILE=-1
set ACTIVEPCT=.1
set ANGLE1=90
set ANGLE2=90
set BOXANG=-45
set NUMPERLAYER=30
set STRESSACTIVATE=0.3
set SCREENSHOT=1
start "SmarticlesSystem_Dynamic.exe" /high "D:\ChronoCode\chronoPkgs\SmarticlesBuild\Release\SmarticlesSystem_Dynamic.exe" %LW% %DT% %NUMLAYER% %READFILE% %ACTIVEPCT% %ANGLE1% %ANGLE2% %BOXANG% %NUMPERLAYER% %STRESSACTIVATE% %SCREENSHOT%