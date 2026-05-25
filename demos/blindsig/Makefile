.PHONY: all clean

LAZER_DIR = ../..
HEXL_BUILD = $(LAZER_DIR)/third_party/hexl-development/build/hexl
# cmake's GNUInstallDirs picks lib or lib64 depending on the distro; try both.
HEXL_LIB = $(firstword $(wildcard $(HEXL_BUILD)/lib64/libhexl.a)\
 $(wildcard $(HEXL_BUILD)/lib/libhexl.a) $(HEXL_BUILD)/lib64/libhexl.a)

LIBS = $(LAZER_DIR)/liblazer.a -lmpfr -lgmp -lm $(HEXL_LIB) -lstdc++ #XXX
CFLAGS = -Wall -Wextra -Wshadow -Wundef -O3 -g

all: blindsig-demo

blindsig-demo: blindsig-demo.c blindsig-demo.h
	$(CC) $(CFLAGS) -I$(LAZER_DIR) -o $@ blindsig-demo.c $(LIBS)

clean:
	rm -f blindsig-demo
