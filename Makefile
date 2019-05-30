NVCC=/usr/local/cuda-9.0/bin/nvcc
CXX ?= g++
GENCODES ?= 61

INCLUDE_DIRS = -I./src
NVCC_FLAGS = -ccbin $(CXX) -std=c++11 -Xcompiler -Wall,-Wextra
NVCC_OPT_FLAGS = -DNDEBUG
NVCC_TEST_FLAGS = -lineinfo
NVCC_DBG_FLAGS = -g -G
NVCC_LIBS = -lstdc++
NVCC_TEST_LIBS = -lgtest

all:
	@echo "Please run 'make check' or 'make bench'."

tests/test-suite: tests/test-suite.cu
	$(NVCC) $(NVCC_TEST_FLAGS) $(NVCC_FLAGS) $(INCLUDE_DIRS) $(NVCC_LIBS) $(NVCC_TEST_LIBS) -o $@ $<
	#$(NVCC) $(NVCC_TEST_FLAGS) $(NVCC_FLAGS) $(GENCODES:%=--gpu-architecture=compute_%) $(GENCODES:%=--gpu-code=sm_%) $(INCLUDE_DIRS) $(NVCC_LIBS) $(NVCC_TEST_LIBS) -o $@ $<

check: tests/test-suite
	@./tests/test-suite

bench/bench: bench/bench.cu
	#$(NVCC) $(NVCC_OPT_FLAGS) $(NVCC_FLAGS) $(INCLUDE_DIRS) $(NVCC_LIBS) -o $@ $<
	$(NVCC) $(NVCC_OPT_FLAGS) $(NVCC_FLAGS) $(GENCODES:%=--gpu-architecture=compute_%) $(GENCODES:%=--gpu-code=sm_%) $(INCLUDE_DIRS) $(NVCC_LIBS) -o $@ $<

bench: bench/bench

paillier: paillier_test.cu
	$(NVCC) $(NVCC_OPT_FLAGS) $(NVCC_FLAGS) $(GENCODES:%=--gpu-architecture=compute_%) $(GENCODES:%=--gpu-code=sm_%) $(INCLUDE_DIRS) $(NVCC_LIBS) -o $@ $<

.PHONY: clean
clean:
	$(RM) tests/test-suite bench/bench
