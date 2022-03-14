@echo off
set TARGET_ARCH=
if "%PROCESSOR_ARCHITECTURE%"=="AMD64" (set X=X64) else (set X=X86)
if /I "%1"=="x86" (set TARGET_ARCH=X86)
if /I "%1"=="x86_amd64" (set TARGET_ARCH=X64)
if /I "%1"=="x86_x64" (set TARGET_ARCH=X64)
if /I "%1"=="amd64" (set TARGET_ARCH=X64)
if /I "%1"=="x64" (set TARGET_ARCH=X64)
if /I "%1"=="amd64_x86" (set TARGET_ARCH=X86)
if /I "%1"=="x64_x86" (set TARGET_ARCH=X86)
if /I "%1"=="" (set TARGET_ARCH=%X%) else (set X=%1%)
if "%TARGET_ARCH%"=="" (goto :end)
set VCVARSALL=""
if exist "%PROGRAMFILES(X86)%" (
for /f "delims=" %%i in ('where /R "%PROGRAMFILES(X86)%" /Q vcvarsall.bat') do set VCVARSALL="%%i"
)
if %VCVARSALL%=="" (for /f "delims=" %%i in ('where /R "%PROGRAMFILES%" vcvarsall.bat') do set VCVARSALL="%%i")
if %VCVARSALL%=="" (goto :end)
call %VCVARSALL% %X%
echo on
nmake MACHINE=%TARGET_ARCH% -f MakefileMSVC.mak
:end
pause