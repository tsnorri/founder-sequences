include ../local.mk
include ../common.mk

OBJECTS		=	bipartite_matcher.o \
				cmdline.o \
				create_segment_texts_task.o \
				generate_context.o \
				greedy_matcher.o \
				join_context.o \
				main.o \
				merge_segments_task.o \
				segment_text.o \
				segmentation_dp_arg.o \
				segmentation_lp_context.o \
				segmentation_sp_context.o \
				update_pbwt_task.o

all: founder_sequences

all-version-suffix: all
	mv founder_sequences `./founder_sequences --version | sed 's/ /-/g'`

clean:
	$(RM) $(OBJECTS) founder_sequences cmdline.c cmdline.h

founder_sequences: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../lib/libbio/src/libbio.a -ldl

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
