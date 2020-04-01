# Redefine for another C++ compiler
CXX = g++
CXXFLAGS = -Wall -O

NFP = nfwprops.o cosbits.o cospar.o phyconst.o labelval.o integ.o newton.o
CSI = cosinfo.o cosbits.o cospar.o phyconst.o labelval.o integ.o
CC_SCE = cosbits.cc cospar.cc labelval.cc phyconst.cc integ.cc \
	 cosinfo.cc nfwprops.cc newton.cc

ALL = cosinfo nfwprops

# Automate dependencies (together with final line) for g++ and GNU make
# - see info for GNU make.
%.d: %.cc
	@$(CXX) -MM $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@: ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

.cc.o:
	$(CXX) -c $(CXXFLAGS) $<

cosinfo: $(CSI)
	$(CXX) $(CXXFLAGS) $(CSI) -o $@

nfwprops: $(NFP)
	$(CXX) $(CXXFLAGS) $(NFP) -o $@

clean:
	rm -f *.o *.d core* $(ALL) *~


-include $(CC_SCE:.cc=.d)
