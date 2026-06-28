@echo off
cls
echo %~dps0 > cowan.cfg
echo .
echo                             Program RCG Mod 11k
echo       computes angular matrix elements, energy levels, and atomic spectra
echo                                  ÜÜÜÜÜÜÜÜÜ
echo                                  Û START Û
echo                                  ßßßßßßßßß
call mark_time.bat "%~dp0rcg11k.exe"
del cowan.cfg > nul
