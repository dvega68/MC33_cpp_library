#--------------------------------------------------------
CPP      = g++
SOURCE   = TestMC33.cpp
OBJ      = $(SOURCE:.cpp=.o)
LFLTK    = -lfltk_gl -lfltk
LX11R    = -lXext -lX11
LIBS     = -lGLU -lGL -ldl $(LFLTK) $(LX11R) -lm -s
CPPINCS  =
BIN      = TestMC33
OPTIM    = -Ofast -m64 -Wall -Wextra -funroll-loops
CPPFLAGS = $(CPPINCS) -std=c++11 $(OPTIM)
LDFLAGS  = $(LIBS)
RM       = rm -f

.PHONY: all all-before all-after clean clean-custom

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<

all:	all-before $(BIN) all-after

$(BIN): $(OBJ)
	#$(SOURCE)
	$(CPP) -o $(BIN) $(OBJ) $(LDFLAGS)

clean: clean-custom
	$(RM) $(OBJ) $(BIN)
#--------------------------------------------------------
