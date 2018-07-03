# Multiple-machine Makefile 

SHELL = /bin/sh  
#.IGNORE: 
 
#Files 
 
#  fortran 
SRCF = 
 
#  c / c++ 
SRCC =  
SRCCPP = widom.cpp geometry.cpp output_widom.cpp
INC = geometry.h output_widom.h
LIB = "-llammps -lfftw3 -lmpi -lz"
# Definitions 
 
ROOT = widomCC
EXE =  $(ROOT)_$@ 
SRC =  $(SRCF) $(SRCC) $(SRCCPP) 
OBJ =  $(SRCF:.f=.o) $(SRCC:.c=.o) $(SRCCPP:.cpp=.o) 
#MACROS = -DDEBUG 
# Help 
 
help:  
	@echo 'This will use the machine-specific makefile you specify' 
	@echo 'Type "make target" where target is one of:' 
	@echo '     newton   (for UTK Linux Cluster)' 
	@echo '     ubuntu'
 
# Targets 
newton darter ubuntu: 
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi 
	@cp -p $(SRC) $(INC)  Obj_$@ 
	@cp Makefile.$@ Obj_$@/Makefile 
	cd Obj_$@; \
	$(MAKE) "OBJ = $(OBJ)" "INC = $(INC)" USRLIB=$(LIB) "MACROS =$(MACROS)" "EXE = ../$(EXE)" ../$(EXE) 
#	$(MAKE) OBJ="$(OBJ)" INC="$(INC)" EXE="$(EXE)" $(EXE) 
#	@if [ -d Obj_$@ ]; then cd Obj$@; rm $(SRC) $(INC) Makefile*; fi 
 
# Cleans 
 
clean_newton: 
	rm -rf Obj_newton *.out
 
clean_darter: 
	rm -rf Obj_darter *.out
 
clean_ubuntu: 
	rm -rf Obj_ubuntu *.out  

