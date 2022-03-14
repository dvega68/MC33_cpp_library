# Project: libMC33++
NAME     = MC33++
!IFNDEF MACHINE
MACHINE  = x64
!ENDIF
LIBDIR   = lib\$(MACHINE)
SRC      = source\libMC33++.cpp
OBJ      = libMC33++.obj
BIN      = $(LIBDIR)\$(NAME).lib
CPPFLAGS = /GL /W3 /Gy /I".\include" /O2 /Zc:inline /fp:fast /Gd /Oi /MD /EHsc /nologo /Ot
LIBFLAGS = /LTCG /MACHINE:$(MACHINE) /NOLOGO
 
all: $(BIN)

$(OBJ): $(SRC)
	$(CPP) $(CPPFLAGS) /c $(SRC)

$(BIN): $(OBJ)
	@if not exist $(LIBDIR) mkdir $(LIBDIR)
	LIB /OUT:$(BIN) $(LIBFLAGS) $(OBJ)

clean:
	-1 del $(OBJ)
	-1 del $(BIN)
