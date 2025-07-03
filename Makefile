# Makefile for TRAKD project
# Compiles all C++ files into separate executable files

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pthread
LDFLAGS = -pthread

# Directories
SRCDIR = src
BINDIR = bin

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# Executable files (without .cpp extension and with bin/ directory prefix)
EXECUTABLES = $(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%,$(SOURCES))

# Default target
all: $(BINDIR) $(EXECUTABLES)

# Create bin directory
$(BINDIR):
	mkdir -p $(BINDIR)

# Rule for compiling each .cpp file into an executable
$(BINDIR)/%: $(SRCDIR)/%.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

# Clean compiled files
clean:
	rm -rf $(BINDIR)

# Rebuild (clean + build)
rebuild: clean all

# Install dependencies (if needed)
install-deps:
	@echo "No special dependencies required for this project"
	@echo "Make sure you have g++ installed with C++17 support"

# Run tests (example)
test: all
	@echo "Running basic tests..."
	@echo "Compilation check completed successfully"
	@ls -la $(BINDIR)/

# Show help
help:
	@echo "Available targets:"
	@echo "  all          - Compile all executable files (default)"
	@echo "  clean        - Remove all compiled files"
	@echo "  rebuild      - Clean and rebuild everything"
	@echo "  install-deps - Show dependency information"
	@echo "  test         - Run simple check"
	@echo "  help         - Show this help"
	@echo ""
	@echo "Executable files will be created in the bin/ directory:"
	@echo "  bin/kmer_analyzer"
	@echo "  bin/locus_bed_generator" 
	@echo "  bin/LocusBedGeneratorDetailed"
	@echo "  bin/distance_analyzer_detailed"
	@echo "  bin/dispersed_repeat_finder"

# Additional flags for debugging
debug: CXXFLAGS += -g -DDEBUG
debug: all

# Declare phony targets
.PHONY: all clean rebuild install-deps test help debug
