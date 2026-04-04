FC = gfortran
CC = gcc
LD = gfortran
AR = ar

FCFLAGS = -fPIC -O2 -std=legacy -fno-automatic -frecord-marker=4 -fall-intrinsics -fallow-argument-mismatch
CCFLAGS = -fPIC -O2

LDFLAGS = -shared -lm

OBJS_WRAPPER = fint_api.o
OBJS_C = Fmtmake5.o fint_capi.o

OBJS = blkdta000.o fardpc.o fetao.o ffalue.o ffnta.o fiab.o fiad.o \
		finta.o fnfb.o fprset.o ftra.o fuxgt.o fvalum.o fycle.o \
		$(OBJS_WRAPPER)  $(OBJS_C)

SODIR = ..
ARDIR = $(SODIR)
SO_TARGET = $(SODIR)/libgint.so
AR_TARGET = $(ARDIR)/libgint.a

.SUFFIXES: .f .c .o

.f.o:
	$(FC) $(FCFLAGS) -c $<

.c.o:
	$(CC) $(CCFLAGS) -c $<

all: $(SO_TARGET) $(AR_TARGET)

$(SO_TARGET): $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^

$(AR_TARGET): $(OBJS)
	$(AR) rcs $@ $^
	cp fint_capi.h $(ARDIR)

clean:
	-rm -rf $(OBJS) $(SO_TARGET) $(AR_TARGET)
