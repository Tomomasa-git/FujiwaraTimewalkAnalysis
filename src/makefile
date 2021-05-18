CC  = gcc
CXX = g++

CFLAGS = -O2

BINDIR = ../bin
LIBDIR = ../lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS)
CXXLIBS  = $(ROOTLIBS)

#TARGET1 = GetPede   Setting
#OBJS1   = GetPede.o Setting.o
#TARGET1 = DrawTQ   Setting
#OBJS1   = DrawTQ.o Setting.o
#TARGET1 = MakeTree   Setting
#OBJS2   = MakeTree.o Setting.o
TARGET1 = Twlk_Ana   Setting
OBJS1   = Twlk_Ana.o Setting.o

all: $(TARGET1)

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o core $(BINDIR)/*
