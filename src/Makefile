include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				founder_sequences.o \
				main.o \
				merge_segments_task.o \
				generate_founder_sequences.o

all: founder_sequences

clean:
	$(RM) $(OBJECTS) founder_sequences cmdline.c cmdline.h

founder_sequences: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../lib/libbio/src/libbio.a

main.c : cmdline.c


test_rmq: test_rmq.o
	$(CXX) -o $@ test_rmq.o $(LDFLAGS)

test_pbwt: test_pbwt.o
	$(CXX) -o $@ test_pbwt.o $(LDFLAGS)
