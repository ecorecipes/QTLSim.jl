# Build QTSim from the unmodified source in src/ (this Makefile is not part of QTSim).
# Needs a C99 compiler and GSL, the GNU Scientific Library, with gsl-config on the PATH
# (Debian and Ubuntu: libgsl-dev; Homebrew: gsl). From the package root:
#
#     make -C reference/QTSim        # builds reference/QTSim/build/QTSim
#
# Floating-point contraction is off, so that the build does the arithmetic the source
# writes on every platform; validate/gsl_locked.jl relies on that to match QTSim bit for
# bit.

CFLAGS ?= -O2
GSL_CFLAGS := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)
SRC := $(wildcard src/*.c)

build/QTSim: $(SRC) $(wildcard src/*.h)
	mkdir -p build
	$(CC) $(CFLAGS) -ffp-contract=off $(GSL_CFLAGS) -o $@ $(SRC) $(GSL_LIBS) -lm

check:
	shasum -a 256 -c SHA256SUMS 2>/dev/null || sha256sum -c SHA256SUMS

clean:
	rm -rf build

.PHONY: check clean
