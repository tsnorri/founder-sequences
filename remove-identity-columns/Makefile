include ../local.mk
include ../common.mk

OBJECTS		=	cmdline.o \
				main.o

all: remove_identity_columns

clean:
	$(RM) $(OBJECTS) remove_identity_columns cmdline.c cmdline.h config.h

remove_identity_columns: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LDFLAGS) ../lib/libbio/src/libbio.a

remove_identity_columns.cc : cmdline.c
cmdline.c : config.h

include ../config.mk
