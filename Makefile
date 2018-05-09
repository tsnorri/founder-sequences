include local.mk
include common.mk

DEPENDENCIES =	lib/lemon/build/lemon/libemon.a \
				lib/libbio/src/libbio.a

ifeq ($(shell uname -s),Linux)
	DEPENDENCIES += lib/swift-corelibs-libdispatch/build/src/libdispatch.a
endif


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
	$(RM) -rf lib/swift-corelibs-libdispatch/build

dependencies: $(DEPENDENCIES)

lib/lemon/build/lemon/libemon.a:
	$(RM) -rf lib/lemon/build && \
	cd lib/lemon && \
	$(MKDIR) build && \
	cd build && \
	CC="$(CC)" \
	CXX="$(CXX)" \
	CFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CFLAGS)" \
	CXXFLAGS="$(SYSTEM_CPPFLAGS) $(SYSTEM_CXXFLAGS)" \
	LDFLAGS="$(SYSTEM_LDFLAGS)" \
	$(CMAKE) ..
	$(MAKE) -C lib/lemon/build VERBOSE=1
	cd lib/lemon/lemon && $(CP) ../build/lemon/config.h ./

lib/libbio/src/libbio.a:
	$(CP) local.mk lib/libbio
	$(MAKE) -C lib/libbio

lib/swift-corelibs-libdispatch/build/src/libdispatch.a:
	$(RM) -rf lib/swift-corelibs-libdispatch/build && \
	cd lib/swift-corelibs-libdispatch && \
	$(MKDIR) build && \
	cd build && \
	$(CMAKE) -G Ninja -DCMAKE_C_COMPILER="$(CC)" -DCMAKE_CXX_COMPILER="$(CXX)" -DBUILD_SHARED_LIBS=OFF .. && \
	$(NINJA)
