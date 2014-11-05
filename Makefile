FC=ifort
FARGS=-standard-semantics -std08
CDIINC=/home/gzhou/Fakeroot/usr/local/include
CDILIB=/home/gzhou/Fakeroot/usr/local/lib
LIBS=cdi

ttest: ttest.o
	$(FC) -o $@ $< -L$(CDILIB) -l$(LIBS)

ttest.o: ttest.f90
	$(FC) $(FARGS) -c $< -I$(CDIINC)

clean:
	rm *.o ttest
