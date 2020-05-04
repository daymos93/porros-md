MAKEFILE = Makefile
exe = porros
fcomp = gfortran #ifort # /opt/intel/compiler70/ia32/bin/ifc  
# Warning: the debugger doesn't get along with the optimization options
# So: not use -O3 WITH -g option
flags =  -O3   
# Remote compilation
OBJS = ziggurat.o porros.o vars.o subroutines.o

.SUFFIXES:            # this deletes the default suffixes 
.SUFFIXES: .f90 .o    # this defines the extensions I want 

.f90.o:  
	$(fcomp) -c $(flags) $< 
        

$(exe):  $(OBJS) Makefile 
	$(fcomp) $(flags) -o $(exe) $(OBJS) 


clean:
	rm ./*.o ./*.mod	


porros.o: porros.f90 ziggurat.o vars.o subroutines.o
vars.o: vars.f90
subroutines.o: subroutines.f90  
