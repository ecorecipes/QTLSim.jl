# QTSimFast hand-over package: builds the original QTSim (original/src, unmodified) and
# QTSimFast (fast/src, a modified version) and tests the second against the first.
# Needs a C99 compiler with C11 atomics, POSIX threads, make, and GSL (the GNU Scientific
# Library) with gsl-config on the PATH (Debian and Ubuntu: libgsl-dev; Homebrew: gsl).
#
#     make               both programs: original/build/QTSim, fast/build/QTSimFast
#     make check         the original's files have their published checksums, and the
#                        files of fast/src said to be unmodified are unmodified
#     make demo          ./demo_parity.sh: a few cases with checksums and times
#     make parity        ./test_parity.sh quick: a small parity matrix (about a minute)
#     make parity-full   the whole matrix (about 1.5 CPU-hours), with sanitizer builds
#     make model-check   the C routines against the recorded output of the Lean model
#     make bench         timings; writes bench_results.md, .csv, bench_runs.csv
#     make docs          formal/render.sh: regenerate the rendered Lean documents and PROOFS.html/.pdf
#     make clean
#
# Floating-point contraction is off so that both builds do the arithmetic the sources
# write. No option that relaxes IEEE arithmetic (-ffast-math, -Ofast) may be added.

CFLAGS ?= -O2
GSL_CFLAGS := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)
OSRC := $(wildcard original/src/*.c)
OHDR := $(wildcard original/src/*.h)
FSRC := $(wildcard fast/src/*.c)
FHDR := $(wildcard fast/src/*.h)
# Each program is compiled from inside its own directory, as src/*.c: with no argument
# both print the name of their main source file as the compiler was given it (__FILE__),
# and this keeps that line the same in the two programs.
OREL := $(patsubst original/%,%,$(OSRC))
FREL := $(patsubst fast/%,%,$(FSRC))
SHARED := GenLib.c PassOptionFile.c RandLib.c
UNMODIFIED := GenLib.c GenLib.h PassOptionFile.c PassOptionFile.h RandLib.c RandLib.h
FB := fast/build

all: original/build/QTSim $(FB)/QTSimFast

original/build/QTSim: $(OSRC) $(OHDR)
	mkdir -p original/build
	cd original && $(CC) $(CFLAGS) -ffp-contract=off $(GSL_CFLAGS) -o build/QTSim $(OREL) $(GSL_LIBS) -lm

$(FB)/QTSimFast: $(FSRC) $(FHDR)
	mkdir -p $(FB)
	cd fast && $(CC) $(CFLAGS) -ffp-contract=off -pthread $(GSL_CFLAGS) -o build/QTSimFast $(FREL) $(GSL_LIBS) -lm

# Variants for the full parity matrix: sanitizers; GSL's own generator, no threads, no
# vector extensions; the original with the sanitizers.
$(FB)/QTSimFast-asan: $(FSRC) $(FHDR)
	mkdir -p $(FB)
	cd fast && $(CC) -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined -ffp-contract=off -pthread $(GSL_CFLAGS) -o build/QTSimFast-asan $(FREL) $(GSL_LIBS) -lm

$(FB)/QTSimFast-tsan: $(FSRC) $(FHDR)
	mkdir -p $(FB)
	cd fast && $(CC) -O1 -g -fsanitize=thread -ffp-contract=off -pthread $(GSL_CFLAGS) -o build/QTSimFast-tsan $(FREL) $(GSL_LIBS) -lm

$(FB)/QTSimFast-plain: $(FSRC) $(FHDR)
	mkdir -p $(FB)
	cd fast && $(CC) $(CFLAGS) -DQTSIM_GSL_RNG -DQTSIM_NO_THREADS -DQTSIM_NO_VECTOR -ffp-contract=off $(GSL_CFLAGS) -o build/QTSimFast-plain $(FREL) $(GSL_LIBS) -lm

$(FB)/QTSim-asan: $(OSRC) $(OHDR)
	mkdir -p $(FB)
	cd original && $(CC) -O1 -g -fsanitize=address,undefined -ffp-contract=off $(GSL_CFLAGS) -o ../fast/build/QTSim-asan $(OREL) $(GSL_LIBS) -lm

variants: $(FB)/QTSimFast-asan $(FB)/QTSimFast-tsan $(FB)/QTSimFast-plain $(FB)/QTSim-asan

check:
	@cd original && (shasum -a 256 -c SHA256SUMS 2>/dev/null || sha256sum -c SHA256SUMS)
	@for f in $(UNMODIFIED); do cmp fast/src/$$f original/src/$$f || exit 1; done; echo "unmodified files of fast/src: OK"

demo: all
	./demo_parity.sh

parity: all
	./test_parity.sh quick

parity-full: all variants
	./test_parity.sh full

bench: all
	./bench.sh

# The two programs' own routines on the instance of formal/ModelCheck.lean.
$(FB)/model_check_orig: formal/model_check.c $(OSRC) $(OHDR)
	mkdir -p $(FB)
	$(CC) $(CFLAGS) -w -ffp-contract=off -DMODEL_ORIG -Ioriginal/src $(GSL_CFLAGS) -o $@ formal/model_check.c $(addprefix original/src/,$(SHARED)) $(GSL_LIBS) -lm

$(FB)/model_check_fast: formal/model_check.c $(FSRC) $(FHDR)
	mkdir -p $(FB)
	$(CC) $(CFLAGS) -w -ffp-contract=off -pthread -DMODEL_FAST -DQTSIM_GSL_RNG -Ifast/src $(GSL_CFLAGS) -o $@ formal/model_check.c $(addprefix fast/src/,$(SHARED)) fast/src/FastRNG.c fast/src/FastThreads.c $(GSL_LIBS) -lm

model-check: $(FB)/model_check_orig $(FB)/model_check_fast
	$(FB)/model_check_orig > $(FB)/model_check_orig.txt
	$(FB)/model_check_fast > $(FB)/model_check_fast.txt
	cmp $(FB)/model_check_orig.txt formal/model_check_expected.txt
	cmp $(FB)/model_check_fast.txt formal/model_check_expected.txt
	@echo "model check: the original's routines, QTSimFast's routines and the Lean model agree (formal/model_check_expected.txt)"

docs:
	formal/render.sh

clean:
	rm -rf original/build fast/build

.PHONY: all variants check demo parity parity-full bench model-check docs clean
