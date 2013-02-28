@echo off

cd ..\..

REM Check and build
SET PATH=%PATH%;C:\Rtools\bin;c:\Rtools\gcc-4.6.3\bin
Rcmd check lmclus
Rcmd build lmclus

REM Cleanup
del lmclus\src\*.o
del lmclus\src\*.dll
del lmclus\src\*.rds
RD /S /Q lmclus.Rcheck

REM Install
Rcmd INSTALL --build -l %R_LIBS% lmclus_1.0.0.tar.gz