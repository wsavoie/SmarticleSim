$EmailFrom = "wsavoielabemail@gmail.com"
$EmailTo = "willsavoie@gmail.com"
$Subject = "sim complete"
$Body = "sim's done brother"
$SMTPServer = "smtp.gmail.com"
$SMTPClient = New-Object Net.Mail.SmtpClient($SmtpServer, 587)
$SMTPClient.EnableSsl = $true
$SMTPClient.Credentials = New-Object System.Net.NetworkCredential("wsavoielabemail@gmail.com", "crablab123");
$SMTPClient.Send($EmailFrom, $EmailTo, $Subject, $Body)