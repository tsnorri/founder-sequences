include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o \
				match_founder_sequences.o

all: match_founder_sequences

clean:
	$(RM) $(OBJECTS) match_founder_sequences cmdline.c cmdline.h

match_founder_sequences: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../lib/libbio/src/libbio.a

main.c : cmdline.c
