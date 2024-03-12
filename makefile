#compiler
FC=gfortran
SRCDIR=lblibs
OBJ= cgamma.o piN.o Deck.o $(SRCDIR)/kinem_utils.o $(SRCDIR)/dirac_utils.o $(SRCDIR)/init_amps.o $(SRCDIR)/regge.o \
$(SRCDIR)/vegas.o frames.o $(SRCDIR)/gauss_int.o $(SRCDIR)/numerics.o pwave.o swave.o dwave.o main.o
FFLAGS= -O3 -g
main: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
vegas.o: $(SRCDIR)/vegas.f
	$(FC) $(FFLAGS) -c $(SRCDIR)/vegas.f
dirac_utils.o: $(SRCDIR)/dirac_utils.f
	$(FC) $(FFLAGS) -c $(SRCDIR)/dirac_utils.f
frames.o: frames.f95
	$(FC) $(FFLAGS) -c frames.f95
cgamma.o: cgamma.f95
	$(FC) $(FFLAGS) -c cgamma.f95
piN.o: piN.f
	$(FC) $(FFLAGS) -c piN.f	
Deck.o: Deck.f
	$(FC) $(FFLAGS) -c Deck.f
pwave.o: pwave.f95
	$(FC) $(FFLAGS) -c pwave.f95
swave.o: swave.f95
	$(FC) $(FFLAGS) -c swave.f95
dwave.o: dwave.f95
	$(FC) $(FFLAGS) -c dwave.f95
main.o:: main.f95
	$(FC) $(FFLAGS) -c main.f95
clean:
	rm -f main *.o $(SRCDIR)/*.o *.mod*
