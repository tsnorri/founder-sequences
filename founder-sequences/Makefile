include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				generate_founder_sequences.o \
				main.o \
				segment_text.o \
				segmentation_dp_arg.o \
				segmentation_lp_context.o \
				segmentation_sp_context.o \
				update_pbwt_task.o

all: founder_sequences

clean:
	$(RM) $(OBJECTS) founder_sequences cmdline.c cmdline.h

founder_sequences: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../lib/libbio/src/libbio.a

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
