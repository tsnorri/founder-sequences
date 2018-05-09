# Build a mostly static binary.

CC					= clang-5.0
CXX					= clang++-5.0
MKDIR				= mkdir
GENGETOPT			= gengetopt

CFLAGS				= -fblocks
CXXFLAGS			= -fblocks

BOOST_ROOT			= /home/tnorri/local/boost-1-67-0-clang++-5.0.2-libstdc++
BOOST_INCLUDE		= -I$(BOOST_ROOT)/include
BOOST_LIBS			= -L$(BOOST_ROOT)/lib -lboost_iostreams
LIBDISPATCH_LIBS	= ../lib/swift-corelibs-libdispatch/libdispatch-build/src/libdispatch.a /usr/lib/x86_64-linux-gnu/libkqueue.a /usr/lib/x86_64-linux-gnu/libBlocksRuntime.a -lbsd -lpthread -static-libstdc++ -static-libgcc
