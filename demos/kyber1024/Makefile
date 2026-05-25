.PHONY: all clean

LAZER_DIR = ../..
HEXL_BUILD = $(LAZER_DIR)/third_party/hexl-development/build/hexl
# cmake's GNUInstallDirs picks lib or lib64 depending on the distro; try both.
HEXL_LIB = $(firstword $(wildcard $(HEXL_BUILD)/lib64/libhexl.a)\
 $(wildcard $(HEXL_BUILD)/lib/libhexl.a) $(HEXL_BUILD)/lib64/libhexl.a)

LIBS = $(LAZER_DIR)/liblazer.a -lmpfr -lgmp -lm $(HEXL_LIB) -lstdc++ #XXX
CFLAGS = -Wall -Wextra -Wshadow -Wundef -O3 -g

all: kyber1024-demo

kyber1024-demo: params.h kyber1024-demo.c
	$(CC) $(CFLAGS) -I$(LAZER_DIR) -o $@ kyber1024-demo.c $(LIBS)

params.h: params.py
	cd ../../scripts && sage lin-codegen.sage ../demos/kyber1024/params.py > ../demos/kyber1024/params.h

clean:
	rm -f kyber1024-demo