#--------------------------------------------------------
CPP      = g++
SOURCE   = TestMC33_glut.cpp
OBJ      = $(SOURCE:.cpp=.o)
LIBS     = -static -lstdc++ -lfreeglut -lglu32 -lopengl32 -lgdi32 -lwinmm -lm -m64 -s
CPPINCS  =
BIN      = TestMC33_glut
OPTIM		 = -Ofast -m64 -Wall -Wextra -funroll-loops
CPPFLAGS = $(CPPINCS) -DFREEGLUT_STATIC $(OPTIM)
RM       = rm -f

.PHONY: all all-before all-after clean clean-custom

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<

all:	all-before $(BIN) all-after

$(BIN): $(OBJ)
	#$(SOURCE)
	$(CPP) -o $(BIN) $(OBJ) $(LIBS)

clean: clean-custom
	$(RM) $(OBJ) $(BIN)
#--------------------------------------------------------
