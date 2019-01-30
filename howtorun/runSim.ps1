param (
    [int]$v = 1, #if not specified saves in v1 folder 
    [int]$s = 1, #if not specified assumes saving screenshot
    [string]$procedureName="SCFrac"
    )
[string]$mainFold=$PSScriptRoot
[string]$runFile = "$mainFold\SmarticlesSystem_Dynamic.exe"
[string]$procedureFile=$procedureName+".txt"
[string]$dataOutFolder="$mainFold\\..\\$procedureName"
[string]$procedure="$mainFold\\$procedureFile"

#first input nums actuationStart,actuation_amp, actuationSpd,smartStr,Bucket type;
#actuationAmp= units of smartWids
#actuationSpd= units of smartWids/sec
#smartStr = smarticle units
#buckettype=KNOBCYLINDER, HOOKRAISE, STRESSSTICK, CYLINDER, BOX, HULL, FLATHOPPER, HOPPER, DRUM, BOXDROP, HOOKRAISE2,HOOKFRACTURE
#cylinder= 3, hookraise2=10, hookFracture=11

[double[]]$lws=@(0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1)
[int[]]$nl=@(38,30,28,27,26,25,24,24)
[int[]]$npl=@(13,13,11,10,10,9,8,7)
[int[]]$vibAmp=@(5)

#[double[]]$lws=@(1)
#[int[]]$nl=@(1)
#[int[]]$npl=@(1)
#[int[]]$vibAmp=@(0)

[double]$dtVal=0.0002

[int]$re=-1
[int]$pa=1
[int]$ang1=90
[int]$ang2=90
[int]$windx=0
[int]$windy=0

echo $dataOutFolder
#make directory for data
New-Item -ItemType directory -Path $dataOutFolder

for($i=0; $i -lt $lws.length ; $i++)
{
    for($j=0; $j -lt $vibAmp.length; $j++)
    {
        $foldOut="$dataOutFolder\\lw_$($lws[$i])_nl_$($nl[$i])_npl_$($npl[$i])_vib_$($vibAmp[$j])"
        New-Item -ItemType directory -Path $foldOut
        $versFold="$foldOut\v$v"
        New-Item -ItemType directory -Path $versFold
        Copy-Item -Path smarticleMoves.csv -Destination $versFold
        Copy-Item -Path $procedure -Destination $versFold
        Set-Location -Path $versFold

        & $runFile $($lws[$i]) $dtVal $($nl[$i]) $re $pa $ang1 $ang2 0 $($npl[$i]) 0.2 $s 0 $windx $windy 1 $($vibAmp[$j]) $procedure
        Set-Location -Path $mainFold
    }
}


& .\sendemail.ps1 $procedureName