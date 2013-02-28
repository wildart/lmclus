@echo off
SET PATH=%PATH%;%MINGW%\bin
REM SET PATH=%PATH%;c:\Rtools\gcc-4.6.3\
REM cd ..
del CMakeCache.txt
cmake -G "MinGW Makefiles" .
mingw32-make.exe clean
mingw32-make.exe
mingw32-make.exe test
