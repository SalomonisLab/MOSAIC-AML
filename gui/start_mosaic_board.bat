@echo off
setlocal
title MOSAIC-AML decision board

REM =====================================================================
REM  Launch the MOSAIC-AML decision board.
REM  It runs on the cluster and opens here through an SSH tunnel.
REM  Just double-click this file. Keep the window open while you use it.
REM =====================================================================

set "SSH_HOST=bmiclusterp-head"
set "PORT=8766"
set "GUI_DIR=/data/salomonis-archive/LabFiles/Nicholas/AML-multimodal/gui"
set "URL=http://localhost:%PORT%/"
set "SSHOPTS=-o StrictHostKeyChecking=accept-new -o ConnectTimeout=10"

echo(
echo   MOSAIC-AML decision board
echo   cluster : %SSH_HOST%      local port : %PORT%
echo   -------------------------------------------------------------
echo(
echo   [1/3] Making sure the board is running on the cluster...
ssh %SSHOPTS% %SSH_HOST% "pgrep -f '[g]ui_server.py %PORT%' >/dev/null 2>&1 && echo already-running || (cd %GUI_DIR% && BROWSER=none nohup /usr/bin/python3 gui_server.py %PORT% >/tmp/matrixgui.log 2>&1 </dev/null & sleep 2 && echo started)"
if errorlevel 1 goto :ssherror

echo   [2/3] Your browser will open at %URL% shortly...
start "" /min powershell -NoProfile -WindowStyle Hidden -Command "Start-Sleep -Seconds 4; Start-Process '%URL%'"

echo   [3/3] Opening the secure tunnel.
echo(
echo   ==================  KEEP THIS WINDOW OPEN  ==================
echo   The board is available while this window stays open.
echo   Close it (or press Ctrl+C) to disconnect.
echo   The cluster server keeps running after you disconnect; to
echo   stop it:  ssh %SSH_HOST% "pkill -f 'gui_server.py %PORT%'"
echo   ============================================================
echo(
ssh -N %SSHOPTS% -o ServerAliveInterval=30 -L %PORT%:127.0.0.1:%PORT% %SSH_HOST%

echo(
echo   Tunnel closed. The board is no longer reachable from this PC.
pause
exit /b 0

:ssherror
echo(
echo   ERROR: could not reach the cluster over SSH (host "%SSH_HOST%").
echo   Test it in a terminal with:   ssh %SSH_HOST%
echo   If that host name is wrong, edit SSH_HOST at the top of this file.
echo(
pause
exit /b 1
