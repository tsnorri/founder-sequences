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

# Use the name config.h since cmdline.c includes it automatically.
config.h: SHELL := /bin/bash
config.h: Makefile ../.git
	IFS=$$'\n'; \
	TAGS=`git tag -l --points-at HEAD`; \
	HASH=`git rev-parse --short --verify HEAD`; \
	VERSION="DEV"; \
	if [ -n "$${TAGS}" ]; \
	then \
		for name in $${TAGS}; \
		do \
			if [[ "$${name}" == v* ]]; \
			then \
				VERSION="$${name:1}"; \
				break; \
			fi; \
		done; \
	fi; \
	printf "#define CMDLINE_PARSER_VERSION \"%s %s\"\n" "$${VERSION}" "$${HASH}" > config.h
