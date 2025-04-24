@ECHO OFF
ECHO RUN PROGRAM INITIALIZED

REM Step 2: Launch the MASTER configuration in a new window
start "Master" "C:\Program Files\bladeRF\x64\bladeRF-cli.exe" -d "*:serial=608d" -s master.txt -i

ECHO OPERATION COMPLETED
