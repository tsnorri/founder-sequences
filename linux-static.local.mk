# Build a mostly static binary.

CC						= clang-6.0
CXX						= clang++-6.0
LIBDISPATCH_CFLAGS		= -U__STDC_HOSTED__ -isystem /usr/lib/llvm-6.0/lib/clang/6.0.0/include
LIBDISPATCH_CXXFLAGS	= -U__STDC_HOSTED__ -isystem /usr/lib/llvm-6.0/lib/clang/6.0.0/include

CFLAGS					= -fblocks -U__STDC_HOSTED__ -isystem /usr/lib/llvm-6.0/lib/clang/6.0.0/include
CXXFLAGS				= -fblocks -U__STDC_HOSTED__ -isystem /usr/lib/llvm-6.0/lib/clang/6.0.0/include
CPPFLAGS				= -DNDEBUG

BOOST_ROOT				= /home/tnorri/local/boost-1-67-0-clang++-5.0.2-libstdc++
BOOST_INCLUDE			= -I$(BOOST_ROOT)/include
BOOST_LIBS				= -L$(BOOST_ROOT)/lib -lboost_iostreams -lboost_serialization -lboost_system
LIBDISPATCH_LIBS		= ../lib/swift-corelibs-libdispatch/build/src/libdispatch.a /usr/lib/x86_64-linux-gnu/libkqueue.a /usr/lib/x86_64-linux-gnu/libBlocksRuntime.a -lbsd -lpthread -static-libstdc++ -static-libgcc
