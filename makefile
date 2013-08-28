# Compiler options
FC 		:= 	gfortran
BINDIR		:=	bin
VPATH		:=	src
FFLAGS		:=	-O3 -g -C -Wall -Warray-bounds -ffree-form -ffixed-line-length-none -fbounds-check -J$(BINDIR) #-I$(BINDIR)

# Information about this run
INFILE 		:= 	input.in

# All modules
OBJS		:= $(BINDIR)/random.o $(BINDIR)/parameters.o $(BINDIR)/fparser.o $(BINDIR)/io.o $(BINDIR)/physics.o $(BINDIR)/eedf.o $(BINDIR)/interpolation.o $(BINDIR)/ratecoeffs.o $(BINDIR)/particles.o $(BINDIR)/populations.o

# Default rule
all: runner | $(BINDIR)

# Set some make specials
.SUFFIXES:
.SUFFIXES: .f .o .mod 

.PHONY: setid getid

# Build rules

$(BINDIR)/%.o: $(VPATH)/%.f | $(BINDIR)
	$(FC) $(FFLAGS) -c $^ -o $@

$(BINDIR)/%.mod: 

runner: $(OBJS)

$(BINDIR):
	@mkdir $(BINDIR)

clean:
	@echo -n "Cleaning..."
	@rm -rf $(BINDIR) *.mod runner* *.dat
	@echo "done!"
