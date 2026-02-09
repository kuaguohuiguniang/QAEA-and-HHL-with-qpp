# Compiler Settings
CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Iinclude

# Executable Names
TARGET_NOZZLE := nozzle_exec
TARGET_HHL    := hhl_exec

# ---------------------------------------------------------
# Source Definitions
# ---------------------------------------------------------

# Sources for the Nozzle Solver
# Main entry: nozzle.cpp (in root)
NOZZLE_SRCS := \
    nozzle.cpp \
    source/nozzle_app.cpp \
    source/quantum_ode_solver.cpp \
    source/interval_stepper.cpp \
    source/quasi1d_driver.cpp \
    source/boundary_conditions.cpp \
    source/knot_policy.cpp \
    source/mean_integral.cpp \
    source/classical_mean_estimator.cpp \
    source/qpp_mean_estimator.cpp \
    source/steady_state_ic.cpp

# Sources for the HHL Solver
# Main entry: HHL.cpp (in root)
HHL_SRCS := \
    HHL.cpp \
    source/HHLSolver.cpp \
    source/HHLinput.cpp

# Generate Object Lists
NOZZLE_OBJS := $(NOZZLE_SRCS:.cpp=.o)
HHL_OBJS    := $(HHL_SRCS:.cpp=.o)

# ---------------------------------------------------------
# Build Rules
# ---------------------------------------------------------

# Default target: compile both if user types just 'make'
all: nozzle HHL

# --- Nozzle Target ---
nozzle: $(TARGET_NOZZLE)

$(TARGET_NOZZLE): $(NOZZLE_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# --- HHL Target ---
HHL: $(TARGET_HHL)

$(TARGET_HHL): $(HHL_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# --- Generic Compilation Rule ---
# Compiles any .cpp to .o (works for both root and source/ files)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# --- Cleanup ---
clean:
	rm -f $(NOZZLE_OBJS) $(HHL_OBJS) $(TARGET_NOZZLE) $(TARGET_HHL)

# Declare phony targets (not actual files)
.PHONY: all clean nozzle HHL
