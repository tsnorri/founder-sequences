include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o

all: insert_identity_columns

clean:
	$(RM) $(OBJECTS) insert_identity_columns cmdline.c cmdline.h config.h

insert_identity_columns: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../lib/libbio/src/libbio.a -ldl

main.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
