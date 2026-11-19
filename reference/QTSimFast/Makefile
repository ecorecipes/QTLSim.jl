# Build QTSimFast, a modified, faster version of QTSim whose output is bit-identical to
# the original's (see README.md). Needs a C99 compiler with C11 atomics, POSIX threads and
# GSL, the GNU Scientific Library, with gsl-config on the PATH. From the package root:
#
#     make -C reference/QTSimFast          # builds reference/QTSimFast/build/QTSimFast
#     make -C reference/QTSimFast check    # the unmodified files are unmodified
#     make -C reference/QTSimFast demo     # a few cases with checksums and times
#     make -C reference/QTSimFast parity   # quick parity test against reference/QTSim
#     make -C reference/QTSimFast parity-full   # the whole parity matrix (hours)
#     make -C reference/QTSimFast bench    # timings, writes bench_results.{md,csv}
#
# Floating-point contraction is off, as in the original's Makefile, so that the build does
# the arithmetic the source writes. No option that relaxes IEEE arithmetic may be added.

CFLAGS ?= -O2
GSL_CFLAGS := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)
SRC := $(wildcard src/*.c)
HDR := $(wildcard src/*.h)
UNMODIFIED := GenLib.c GenLib.h PassOptionFile.c PassOptionFile.h RandLib.c RandLib.h

all: build/QTSimFast

build/QTSimFast: $(SRC) $(HDR)
	mkdir -p build
	$(CC) $(CFLAGS) -ffp-contract=off -pthread $(GSL_CFLAGS) -o $@ $(SRC) $(GSL_LIBS) -lm

# Variants for testing: sanitizers, GSL's own generator, no threads, no vector extensions.
build/QTSimFast-asan: $(SRC) $(HDR)
	mkdir -p build
	$(CC) -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined -ffp-contract=off -pthread $(GSL_CFLAGS) -o $@ $(SRC) $(GSL_LIBS) -lm

build/QTSimFast-tsan: $(SRC) $(HDR)
	mkdir -p build
	$(CC) -O1 -g -fsanitize=thread -ffp-contract=off -pthread $(GSL_CFLAGS) -o $@ $(SRC) $(GSL_LIBS) -lm

build/QTSimFast-plain: $(SRC) $(HDR)
	mkdir -p build
	$(CC) $(CFLAGS) -DQTSIM_GSL_RNG -DQTSIM_NO_THREADS -DQTSIM_NO_VECTOR -ffp-contract=off $(GSL_CFLAGS) -o $@ $(SRC) $(GSL_LIBS) -lm

# The original, built with the sanitizers, to look for undefined behaviour in it.
build/QTSim-asan: $(wildcard ../QTSim/src/*.c) $(wildcard ../QTSim/src/*.h)
	mkdir -p build
	$(CC) -O1 -g -fsanitize=address,undefined -ffp-contract=off $(GSL_CFLAGS) -o $@ $(wildcard ../QTSim/src/*.c) $(GSL_LIBS) -lm

variants: build/QTSimFast-asan build/QTSimFast-tsan build/QTSimFast-plain build/QTSim-asan

../QTSim/build/QTSim:
	$(MAKE) -C ../QTSim

# The files that were not modified are byte-identical to the original's.
check:
	@for f in $(UNMODIFIED); do cmp src/$$f ../QTSim/src/$$f || exit 1; done; echo "unmodified files: OK"

# A short demonstration with checksums (writes nothing; see PARITY_DEMO.txt for a transcript).
demo: build/QTSimFast ../QTSim/build/QTSim
	./demo_parity.sh

parity: build/QTSimFast ../QTSim/build/QTSim
	./test_parity.sh quick

parity-full: build/QTSimFast ../QTSim/build/QTSim variants
	./test_parity.sh full

bench: build/QTSimFast ../QTSim/build/QTSim
	./bench.sh

# The two programs' own routines on the instance of formal/ModelCheck.lean (see PROOFS.md).
SHARED := GenLib.c PassOptionFile.c RandLib.c
build/model_check_orig: formal/model_check.c $(wildcard ../QTSim/src/*.c) $(wildcard ../QTSim/src/*.h)
	mkdir -p build
	$(CC) $(CFLAGS) -w -ffp-contract=off -DMODEL_ORIG -I../QTSim/src $(GSL_CFLAGS) -o $@ formal/model_check.c $(addprefix ../QTSim/src/,$(SHARED)) $(GSL_LIBS) -lm

build/model_check_fast: formal/model_check.c $(SRC) $(HDR)
	mkdir -p build
	$(CC) $(CFLAGS) -w -ffp-contract=off -pthread -DMODEL_FAST -DQTSIM_GSL_RNG -Isrc $(GSL_CFLAGS) -o $@ formal/model_check.c $(addprefix src/,$(SHARED)) src/FastRNG.c src/FastThreads.c $(GSL_LIBS) -lm

model-check: build/model_check_orig build/model_check_fast
	build/model_check_orig > build/model_check_orig.txt
	build/model_check_fast > build/model_check_fast.txt
	cmp build/model_check_orig.txt formal/model_check_expected.txt
	cmp build/model_check_fast.txt formal/model_check_expected.txt
	@echo "model check: the original's routines, QTSimFast's routines and the Lean model agree"

# The hand-over package (dist/ beside the repository; see make_package.sh).
package:
	./make_package.sh

clean:
	rm -rf build

.PHONY: all check demo model-check package parity parity-full bench variants clean
