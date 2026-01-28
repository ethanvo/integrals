# Makefile for Quantum Chemistry Integral Engine

CC = gcc
CFLAGS = -Wall -Wextra -O2 -I./include
LDFLAGS = -lm

SRC_DIR = src
TEST_DIR = tests
EXAMPLES_DIR = examples
OBJ_DIR = obj
BIN_DIR = bin

# Source files
SOURCES = $(wildcard $(SRC_DIR)/*.c)
OBJECTS = $(SOURCES:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

# Test files
TEST_SOURCES = $(wildcard $(TEST_DIR)/*.c)
TEST_BINS = $(TEST_SOURCES:$(TEST_DIR)/%.c=$(BIN_DIR)/%)

# Example files
EXAMPLE_SOURCES = $(wildcard $(EXAMPLES_DIR)/*.c)
EXAMPLE_BINS = $(EXAMPLE_SOURCES:$(EXAMPLES_DIR)/%.c=$(BIN_DIR)/%)

# Library
LIBRARY = libintegrals.a

.PHONY: all clean tests examples

all: $(OBJ_DIR) $(BIN_DIR) $(LIBRARY)

# Create directories
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Build object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Build static library
$(LIBRARY): $(OBJECTS)
	ar rcs $@ $^

# Build tests
tests: $(TEST_BINS)

$(BIN_DIR)/%: $(TEST_DIR)/%.c $(LIBRARY)
	$(CC) $(CFLAGS) $< -L. -lintegrals $(LDFLAGS) -o $@

# Build examples
examples: $(EXAMPLE_BINS)

$(BIN_DIR)/%: $(EXAMPLES_DIR)/%.c $(LIBRARY)
	$(CC) $(CFLAGS) $< -L. -lintegrals $(LDFLAGS) -o $@

# Clean
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) $(LIBRARY)

# Run tests
run-tests: tests
	@for test in $(TEST_BINS); do \
		echo "Running $$test..."; \
		$$test || exit 1; \
	done
	@echo "All tests passed!"
