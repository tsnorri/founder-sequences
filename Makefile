include local.mk
include common.mk

DEPENDENCIES =	lib/lemon/build/lemon/libemon.a \
				lib/libbio/src/libbio.a


.PHONY: all clean-all clean clean-dependencies dependencies

all: dependencies
	$(MAKE) -C founder-sequences all
	$(MAKE) -C match-sequences-to-founders all

clean-all: clean clean-dependencies

clean:
	$(MAKE) -C founder-sequences clean
	$(MAKE) -C match-sequences-to-founders clean

clean-dependencies:
	$(RM) -rf lib/lemon/build
	$(MAKE) -C lib/libbio clean-all

dependencies: $(DEPENDENCIES)

lib/lemon/build/lemon/libemon.a:
	rm -rf lib/lemon/build && \
	cd lib/lemon && \
	mkdir build && \
	cd build && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	CFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CFLAGS)" \
	CXXFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CXXFLAGS)" \
	LDFLAGS="$(SYSTEM_LDFLAGS)" \
	cmake ..
	$(MAKE) -C lib/lemon/build VERBOSE=1
	cd lib/lemon/lemon && cp ../build/lemon/config.h ./

lib/libbio/src/libbio.a:
	$(CP) local.mk lib/libbio
	$(MAKE) -C lib/libbio
