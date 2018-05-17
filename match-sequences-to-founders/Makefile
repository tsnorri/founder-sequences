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
