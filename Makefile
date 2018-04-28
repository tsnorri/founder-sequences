include local.mk
include common.mk

DEPENDENCIES = lib/lemon/build/lemon/libemon.a


.PHONY: all clean-all clean clean-dependencies dependencies

all: dependencies
	$(MAKE) -C src

clean-all: clean clean-dependencies

clean:
	$(MAKE) -C src clean

clean-dependencies:
	$(RM) -r lib/lemon/build/lemon/libemon.a

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
