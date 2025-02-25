@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2017a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2017a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=fastode45_MEE_mex
set MEX_NAME=fastode45_MEE_mex
set MEX_EXT=.mexw64
call "C:\PROGRA~1\MATLAB\R2017a\sys\lcc64\lcc64\mex\lcc64opts.bat"
echo # Make settings for fastode45_MEE > fastode45_MEE_mex.mki
echo COMPILER=%COMPILER%>> fastode45_MEE_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> fastode45_MEE_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> fastode45_MEE_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> fastode45_MEE_mex.mki
echo LINKER=%LINKER%>> fastode45_MEE_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> fastode45_MEE_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> fastode45_MEE_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> fastode45_MEE_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> fastode45_MEE_mex.mki
echo BORLAND=%BORLAND%>> fastode45_MEE_mex.mki
echo OMPFLAGS= >> fastode45_MEE_mex.mki
echo OMPLINKFLAGS= >> fastode45_MEE_mex.mki
echo EMC_COMPILER=lcc64>> fastode45_MEE_mex.mki
echo EMC_CONFIG=optim>> fastode45_MEE_mex.mki
"C:\Program Files\MATLAB\R2017a\bin\win64\gmake" -B -f fastode45_MEE_mex.mk
