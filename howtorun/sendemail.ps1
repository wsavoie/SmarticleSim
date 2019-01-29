param (
    [string]$simType = "" #if not specified saves in v1 folder 
    )
$EmailFrom = "wsavoielabemail@gmail.com"
$EmailTo = "wsavoie@gatech.edu"
$Subject = "sim $simType complete on $Env:COMPUTERNAME"
$Body = "running sim is complete!"
$SMTPServer = "smtp.gmail.com"
$SMTPClient = New-Object Net.Mail.SmtpClient($SmtpServer, 587)
$SMTPClient.EnableSsl = $true
$SMTPClient.Credentials = New-Object System.Net.NetworkCredential("wsavoielabemail@gmail.com", "crablab123");
$SMTPClient.Send($EmailFrom, $EmailTo, $Subject, $Body)