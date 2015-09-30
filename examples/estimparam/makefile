SRCS=$(wildcard *.cpp)
BINS=$(SRCS:.cpp=)

CXXFLAGS := $(shell pkg-config --define-variable=home=$(HOME) --cflags ibex) 
LIBS	 := $(shell pkg-config --define-variable=home=$(HOME) --libs  ibex)
LIBDIR	 := $(shell pkg-config --define-variable=home=$(HOME) --libdir  ibex) 

ifeq ($(DEBUG), yes)
CXXFLAGS := $(CXXFLAGS) -O0 -g -Wall -frounding-math  -m32
#CXXFLAGS := $(CXXFLAGS) -O0 -g -Wall -frounding-math  
else
#CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG -Wno-deprecated -frounding-math -m32 -pg
CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG -Wno-deprecated -frounding-math -m32 
#CXXFLAGS    := $(CXXFLAGS) -O3 -DNDEBUG -Wno-deprecated -frounding-math 
endif

all: $(BINS)

% :	%.cpp
	echo $(BINS) $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f $(BINS)

