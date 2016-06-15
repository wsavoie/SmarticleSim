set LW=0.8
set DT=.00025
set NUMLAYER=1
REM 2=R&W 1=R 0='' -1=W
set READFILE=0
set ACTIVEPCT=1
set ANGLE1=0
set ANGLE2=0
set BOXANG=-25
set NUMPERLAYER=1
set LAZY=.00025
start "SmarticlesSystem_Dynamic.exe" /high "D:\ChronoCode\chronoPkgs\SmarticlesBuild\Release\SmarticlesSystem_Dynamic.exe" %LW% %DT% %NUMLAYER% %READFILE% %ACTIVEPCT% %ANGLE1% %ANGLE2% %BOXANG% %NUMPERLAYER% %LAZY%