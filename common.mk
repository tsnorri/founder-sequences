# Ignore Xcode's setting since the SDK may contain older versions of Clang and libc++.
unexport SDKROOT

WARNING_FLAGS	= -Wall -Werror -Wno-deprecated-declarations -Wno-unused
OPT_FLAGS		= -O2 -g
#OPT_FLAGS		= -O0 -g

CP				?= cp
MKDIR			?= mkdir
CMAKE			?= cmake
NINJA			?= ninja
GENGETOPT		?= gengetopt
RAGEL			?= ragel
DOT				?= dot
TAR				?= tar

CFLAGS			?=
CXXFLAGS		?=
CPPFLAGS		?=
LDFLAGS			?=
SYSTEM_CFLAGS	?=
SYSTEM_CXXFLAGS	?=
SYSTEM_CPPFLAGS	?=
SYSTEM_LDFLAGS	?=

CFLAGS			+= -std=c99   $(OPT_FLAGS) $(WARNING_FLAGS) $(SYSTEM_CFLAGS)
CXXFLAGS		+= -std=c++17 $(OPT_FLAGS) $(WARNING_FLAGS) $(SYSTEM_CXXFLAGS)
CPPFLAGS		+= $(SYSTEM_CPPFLAGS) -I../include -I../lib/lemon -I../lib/libbio/include -I../lib/sdsl-lite/include -I../lib/swift-corelibs-libdispatch $(BOOST_INCLUDE)
LDFLAGS			+= $(SYSTEM_LDFLAGS) $(BOOST_LIBS) ../lib/lemon/build/lemon/libemon.a

ifeq ($(shell uname -s),Linux)
	LDFLAGS		+= ../lib/swift-corelibs-libdispatch/build/src/libdispatch.a ../lib/swift-corelibs-libdispatch/build/libBlocksRuntime.a -lbsd -lpthread -lz
endif

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

%.c: %.ggo
	$(GENGETOPT) --input="$<"

%.cc: %.rl
	$(RAGEL) -L -C -G2 -o $@ $<

%.dot: %.rl
	$(RAGEL) -V -p -o $@ $<

%.pdf: %.dot
	$(DOT) -Tpdf $< > $@
