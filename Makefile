#################################################################
# Makefile created using the tool 'Creamake'
# 
# Creamake is distributed under the GNU GPL license
# Author: Francisco Pena, fran(dot)pena(at)usc(dot)es
# Download page: http://sourceforge.net/projects/creamake/
#################################################################
 
#################################################################
# User-modifiable options
#################################################################
 
# SOURCE AND COMMONS FOLDERS (separated by spaces)
dir_fuentes = src
 
# OBJECT AND .MOD FOLDER
dir_objetos = object
 
# MAIN SOURCE FILE (include relative path from folder where Makefile is)
condir_principal =         src/main.f90
 
# EXECUTABLE NAME 
ejecutable = firm_finance.exe
 
# NEEDED TO convert ejecutable THE DEFAULT RULE: 
$(ejecutable): $(condir_principal) 
 
# MODULES
modulos = global_variables.f90 functions.f90 LinInterpModule.f90 \
tauchen_mod.f90 grids.f90 solve_values.f90 policy_to_mm.f90 \
calibration_estimation.f90
 
# MODULE DEPENDENCIES
# if pru1 depends on pru2... pru1.o: pru2.o
functions.o: global_variables.o
grids.o: global_variables.o tauchen_mod.o LinInterpModule.o functions.o
solve_values.o: global_variables.o functions.o grids.o LinInterpModule.o
policy_to_mm.o: global_variables.o LinInterpModule.o functions.o \
solve_values.o tauchen_mod.o grids.o
calibration_estimation.o: global_variables.o functions.o grids.o \
LinInterpModule.o solve_values.o policy_to_mm.o
 
# INCLUDES
includes = 
 
# COMPILER
FC = mpif90
 
# COMPILER OPTIONS
FFLAGS = -O3
 
# LINKER OPTIONS
LDFLAGS = -O3
 
#################################################################
# Non-modifiable part
#################################################################
 
# SOURCE FOLDERS
VPATH =   $(subst ,:,$(strip $(dir_fuentes)))
vpath %.o $(dir_objetos)
 
# SOURCES
fuentes_ = $(filter %.f %.F %.for %.FOR %.f90 %.F90 %.f95 %.F95 %.f03 %.F03,$(shell ls $(dir_fuentes)))
fuentes  = $(filter-out $(notdir $(condir_principal)) $(modulos),$(fuentes_))
 
# OBJECTS
modulos_obj = $(addsuffix .o,$(basename $(modulos)))
fuentes_obj = $(addsuffix .o,$(basename $(fuentes)))
 
# OBJECTS WITH PATH
condir_modulos_obj = $(addprefix $(dir_objetos)/,$(modulos_obj))
condir_fuentes_obj = $(addprefix $(dir_objetos)/,$(fuentes_obj))
 
# COMPILATION OPTIONS
FFLAGS += $(patsubst %,-I%,$(dir_fuentes))
FFLAGS += -I$(dir_objetos)
 
# MAIN RULE
all: $(ejecutable)
 
$(ejecutable): $(includes) $(modulos_obj) $(fuentes_obj)
	$(FC) -o $(ejecutable) $(FFLAGS) $(condir_principal) $(condir_modulos_obj) $(condir_fuentes_obj) $(LDFLAGS)
 
# SOURCES RULE
$(fuentes_obj): $(includes) $(modulos_obj)
 
# RULE PATTERNS
%.o:%.f
	$(FC) -c -o $@ $(FFLAGS) $<
	@mv $@ $(dir_objetos) 
%.o:%.F
	$(FC) -c -o $@ $(FFLAGS) $< 
	@mv $@ $(dir_objetos) 
%.o:%.for
	$(FC) -c -o $@ $(FFLAGS) $< 
	@mv $@ $(dir_objetos) 
%.o:%.FOR
	$(FC) -c -o $@ $(FFLAGS) $< 
	@mv $@ $(dir_objetos) 
%.o:%.f90
	$(FC) -c -o $@ $(FFLAGS) $< 
	@mv $@ $(dir_objetos) 
%.o:%.F90
	$(FC) -c -o $@ $(FFLAGS) $< 
	@mv $@ $(dir_objetos) 
 
.PHONY: clean
clean:
	rm $(dir_objetos)/*.o      
	rm $(dir_objetos)/*.mod    
	rm $(ejecutable)
