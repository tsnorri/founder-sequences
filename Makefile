include local.mk
include common.mk

DEPENDENCIES =	lib/lemon/build/lemon/libemon.a \
				lib/libbio/src/libbio.a

ifeq ($(shell uname -s),Linux)
	DEPENDENCIES += lib/swift-corelibs-libdispatch/build/src/libdispatch.a
endif


.PHONY: all clean-all clean clean-dependencies dependencies dist

all: dependencies
	$(MAKE) -C founder-sequences all
	$(MAKE) -C remove-identity-columns all
	$(MAKE) -C insert-identity-columns all
	$(MAKE) -C match-sequences-to-founders all

clean-all: clean clean-dependencies clean-dist

clean:
	$(MAKE) -C founder-sequences clean
	$(MAKE) -C remove-identity-columns clean
	$(MAKE) -C insert-identity-columns clean
	$(MAKE) -C match-sequences-to-founders clean

clean-dependencies: lib/libbio/local.mk
	$(RM) -rf lib/lemon/build
	$(MAKE) -C lib/libbio clean-all
	$(RM) -rf lib/swift-corelibs-libdispatch/build

clean-dist:
	$(RM) -rf founder-sequences-0.1

dist: all
	$(MKDIR) -p founder-sequences-0.1
	$(CP) founder-sequences/founder_sequences founder-sequences-0.1/
	$(CP) remove-identity-columns/remove_identity_columns founder-sequences-0.1/
	$(CP) match-sequences-to-founders/match_founder_sequences founder-sequences-0.1/
	$(CP) README.md founder-sequences-0.1/
	$(CP) LICENSE founder-sequences-0.1/
	$(CP) lib/swift-corelibs-libdispatch/LICENSE founder-sequences-0.1/swift-corelibs-libdispatch-license.txt
	$(TAR) czf founder-sequences-0.1.tar.gz founder-sequences-0.1
	$(RM) -rf founder-sequences-0.1

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

lib/libbio/local.mk:
	$(CP) local.mk lib/libbio

lib/libbio/src/libbio.a: lib/libbio/local.mk
	$(MAKE) -C lib/libbio

lib/swift-corelibs-libdispatch/build/src/libdispatch.a:
	$(RM) -rf lib/swift-corelibs-libdispatch/build && \
	cd lib/swift-corelibs-libdispatch && \
	$(MKDIR) build && \
	cd build && \
	$(CMAKE) \
		-G Ninja \
		-DCMAKE_C_COMPILER="$(CC)" \
		-DCMAKE_CXX_COMPILER="$(CXX)" \
		-DCMAKE_C_FLAGS="$(LIBDISPATCH_CFLAGS)" \
		-DCMAKE_CXX_FLAGS="$(LIBDISPATCH_CXXFLAGS)" \
		-DBUILD_SHARED_LIBS=OFF \
		.. && \
	$(NINJA) -v
