#####################################
#                                   #
#  HADDOCK tools makefile           #
#                                   #
#####################################

CPP=g++
CPFLAGS=-O2

EXEC=contact-segid contact-chainID 

all: 
	make $(EXEC)

contact-segid: contact-segid.cpp
	$(CPP) $(CPFLAGS) -o contact-segid contact-segid.cpp

contact-chainID: contact-chainID.cpp
	$(CPP) $(CPFLAGS) -o contact-chainID contact-chainID.cpp

clean :
	/bin/rm $(EXEC)

