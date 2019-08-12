include local.mk
include common.mk

DEPENDENCIES =	lib/lemon/build/lemon/libemon.a \
				lib/libbio/src/libbio.a
ifeq ($(shell uname -s),Linux)
	DEPENDENCIES += lib/swift-corelibs-libdispatch/build/src/libdispatch.a
endif

# “$() $()” is a literal space.
OS_NAME = $(shell tools/os_name.sh)
VERSION = $(subst $() $(),-,$(shell tools/git_version.sh))
DIST_TARGET_DIR = founder-sequences-$(VERSION)
DIST_NAME_SUFFIX = $(if $(TARGET_TYPE),-$(TARGET_TYPE),)
DIST_TAR_GZ = founder-sequences-$(VERSION)-$(OS_NAME)$(DIST_NAME_SUFFIX).tar.gz


.PHONY: all clean-all clean clean-dependencies dependencies dist

all: $(DEPENDENCIES)
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

dist: $(DIST_TAR_GZ)

$(DIST_TAR_GZ):	founder-sequences/founder_sequences \
				insert-identity-columns/insert_identity_columns \
				match-sequences-to-founders/match_founder_sequences \
				remove-identity-columns/remove_identity_columns
	$(MKDIR) -p $(DIST_TARGET_DIR)
	$(CP) founder-sequences/founder_sequences $(DIST_TARGET_DIR)
	$(CP) insert-identity-columns/insert_identity_columns $(DIST_TARGET_DIR)
	$(CP) match-sequences-to-founders/match_founder_sequences $(DIST_TARGET_DIR)
	$(CP) remove-identity-columns/remove_identity_columns $(DIST_TARGET_DIR)
	$(CP) README.md $(DIST_TARGET_DIR)
	$(CP) LICENSE $(DIST_TARGET_DIR)
	$(CP) lib/swift-corelibs-libdispatch/LICENSE $(DIST_TARGET_DIR)/swift-corelibs-libdispatch-license.txt
	$(TAR) czf $(DIST_TAR_GZ) $(DIST_TARGET_DIR)
	$(RM) -rf $(DIST_TARGET_DIR)

dependencies: $(DEPENDENCIES)

lib/lemon-1.3.1.tar.gz:
	cd lib && $(WGET) http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz

lib/lemon: lib/lemon-1.3.1.tar.gz
	cd lib && $(RM) -rf lemon && $(MKDIR) -p lemon && $(TAR) -xzf lemon-1.3.1.tar.gz -C lemon --strip-components 1

lib/lemon/build/lemon/libemon.a: lib/lemon
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

lib/libbio/local.mk: local.mk
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
