.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
objdir = libtracer/include
bindir = libtracer/bin
libdir = libs
srcdir = src

# Command-line options at make call
env   ?= None      # options: manto,eolo,airaki,iplex
debug ?= 0 

ifeq ($(env),manto) ## env=manto

    ## IFORT OPTIONS ##
    FC  = ifort
    INC_NC  = -I/home/jalvarez/work/librairies/netcdflib/include
    LIB_NC  = -L/home/jalvarez/work/librairies/netcdflib/lib -lnetcdf
    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

    FLAGS    = -module $(objdir) -L$(objdir) $(INC_NC)
    LFLAGS   = $(LIB_NC) $(LIB_MKL)

    DFLAGS   = -vec-report0 -O2 -fp-model precise -i_dynamic 
    ifeq ($(debug), 1)
        DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0 -fp-model precise -i_dynamic 
    endif

else ifeq ($(env),eolo) ## env=eolo

#    ## IFORT OPTIONS ##
#    FC  = ifort
#    INC_NC  = -I/home/fispalma22/work/librairies/netcdflib/include
#    LIB_NC  = -L/home/fispalma22/work/librairies/netcdflib/lib -lnetcdf
#    INC_COORD = -I/home/fispalma25/robinson/models/EURICE/coord/.obj
#    LIB_COORD = /home/fispalma25/robinson/models/EURICE/coord/libcoordinates.a
#    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#
#    FLAGS    = -module $(objdir) -L$(objdir) $(INC_COORD) $(INC_NC)
#    LFLAGS   = $(LIB_COORD) $(LIB_NC) $(LIB_MKL)
#
#    DFLAGS   = -vec-report0 -O2 -fp-model precise
#    ifeq ($(debug), 1)
#        DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0 -fp-model precise
#    endif

    ## GFORTRAN OPTIONS ##
    FC  = gfortran
    INC_NC  = -I/home/fispalma25/apps/netcdf/netcdf/include
    LIB_NC  = -L/home/fispalma25/apps/netcdf/netcdf/lib -lnetcdff -lnetcdf
    INC_COORD = -I/home/fispalma25/apps/coordinates/.obj
    LIB_COORD = /home/fispalma25/apps/coordinates/libcoordinates.a

    FLAGS  = -I$(objdir) -J$(objdir) $(INC_COORD) $(INC_NC)
    LFLAGS = $(LIB_COORD) $(LIB_NC)

    DFLAGS = -O3
    ifeq ($(debug), 1)  # ,underflow
        DFLAGS   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all
    endif

else ifeq ($(env),airaki) ## env=airaki

    ## GFORTRAN OPTIONS ##
    FC  = gfortran
    INC_NC  = -I/opt/local/include
    LIB_NC  = -L/opt/local/lib -lnetcdff -lnetcdf
    INC_COORD = -I/Users/robinson/models/EURICE/coordinates/.obj
    LIB_COORD = /Users/robinson/models/EURICE/coordinates/libcoordinates.a

    FLAGS  = -I$(objdir) -J$(objdir) $(INC_COORD) $(INC_NC) 
    LFLAGS = $(LIB_COORD) $(LIB_NC)

    DFLAGS = -O3
    ifeq ($(debug), 1)  # ,underflow
        DFLAGS   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all
    endif

else ifeq ($(env),pik) ## env=pik

    ## IFORT OPTIONS ##
    FC  = ifort
    INC_NC  = -I${NETCDF_FORTRANROOT}/include
    LIB_NC  = -L${NETCDF_FORTRANROOT}/lib -lnetcdff -L${NETCDF_CROOT}/lib -lnetcdf
    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
    INC_COORD = -I/p/projects/tumble/robinson/EURICE/coord/.obj
	LIB_COORD = /p/projects/tumble/robinson/EURICE/coord/libcoordinates.a

    FLAGS    = -module $(objdir) -L$(objdir) $(INC_COORD) $(INC_NC) 
    LFLAGS   = $(LIB_COORD) $(LIB_NC)

    DFLAGS   = -O3
    ifeq ($(debug), 1)
        DFLAGS   = -C -g -traceback -ftrapuv -fpe0 -check all
    endif

else 
    
    ## None ##
    FC = $(error "Define env")

endif

## Individual libraries or modules ##
$(objdir)/nml.o: $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/ncio.o: $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/tracer_precision.o: $(srcdir)/tracer_precision.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/tracer_interp.o: $(srcdir)/tracer_interp.f90 $(objdir)/tracer_precision.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/tracer3D.o: $(srcdir)/tracer3D.f90 $(objdir)/tracer_precision.o $(objdir)/nml.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/tracer2D.o: $(srcdir)/tracer2D.f90 $(objdir)/tracer3D.o $(objdir)/tracer_precision.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/tracer_io.o: $(srcdir)/tracer_io.f90 $(objdir)/tracer3D.o $(objdir)/tracer_precision.o $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/tracer.o: $(srcdir)/tracer.f90 $(objdir)/tracer_precision.o $(objdir)/nml.o $(objdir)/ncio.o \
					$(objdir)/tracer_interp.o $(objdir)/tracer3D.o $(objdir)/tracer2D.o $(objdir)/tracer_io.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

## bspline-fortran #####
$(objdir)/bspline_sub_module.o: $(libdir)/bspline-fortran/bspline_sub_module.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/bspline_oo_module.o: $(libdir)/bspline-fortran/bspline_oo_module.f90 $(objdir)/bspline_sub_module.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/bspline_module.o: $(libdir)/bspline-fortran/bspline_module.f90 $(objdir)/bspline_sub_module.o $(objdir)/bspline_oo_module.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

obj_bspline =   $(objdir)/bspline_sub_module.o    \
				$(objdir)/bspline_oo_module.o   \
				$(objdir)/bspline_module.o

#########################

obj_tracer =    $(objdir)/tracer_precision.o \
				$(objdir)/tracer_interp.o \
				$(objdir)/tracer3D.o \
				$(objdir)/tracer2D.o \
				$(objdir)/tracer_io.o \
				$(objdir)/tracer.o
				

obj_libs   =    $(objdir)/nml.o    \
				$(objdir)/ncio.o
				
# Static library
tracer-static: $(obj_libs) $(obj_bspline) $(obj_tracer)
	ar rc $(objdir)/libtracer.a $(obj_bspline) $(obj_tracer)
	@echo " "
	@echo "    $(objdir)/libtracer.a is ready."
	@echo " "
			   
## Complete programs

test: tracer-static 
	$(FC) $(DFLAGS) $(FLAGS) -o $(bindir)/test_greenland.x $(srcdir)/test_greenland.f90 \
			-L${CURDIR}/libtracer/include -ltracer $(obj_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_greenland.x is ready."
	@echo " "

test_profile: tracer-static
	$(FC) $(DFLAGS) $(FLAGS) -o $(bindir)/test_profile.x $(srcdir)/test_profile.f90 \
			-L${CURDIR}/libtracer/include -ltracer $(obj_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_profile.x is ready."
	@echo " "

clean:
	rm -f *.x $(bindir)/*.x
	rm -f gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
