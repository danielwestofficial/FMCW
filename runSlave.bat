@ECHO OFF
ECHO RUN PROGRAM INITIALIZED

REM Launch the slave configuration in a new window
start "Slave" "C:\Program Files\bladeRF\x64\bladeRF-cli.exe" -d "*:serial=7789" -s slave.txt -i

ECHO Both Master and Slave have been closed
