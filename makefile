# --- Project config ----------------------------------------------------------
APP           ?= AMDAT
SRC_DIR       := src
BUILD_DIR     := build

# Feature toggles
MODE          ?= release   # debug | release
OMP           ?= 1         # 1 enables OpenMP

# Normalize MODE and OMP to avoid trailing-space bugs
MODE := $(strip $(MODE))
OMP  := $(strip $(OMP))

# Toolchains
CXX           ?= g++
CC            ?= gcc
STD           ?= c++17
WARN          := -Wall -Wextra -Wshadow -Wconversion -Wpedantic

# Wavevector paths
WV3D ?= $(CURDIR)/src/qvectors/qvectors3d/qvector
WV2D ?= $(CURDIR)/src/qvectors/qvectors2d/qvector
WV1D ?= $(CURDIR)/src/qvectors/qvectors1d/qvector

# --- Flags -------------------------------------------------------------------
CPPFLAGS := -MMD -MP -I./src \
						-I./src/voro++-0.4.6/src \
						-L./src/voro++-0.4.6/src \
						-I./xdrfile-1.1b/src \
            -DWV3D=\"$(WV3D)\" \
            -DWV2D=\"$(WV2D)\" \
            -DWV1D=\"$(WV1D)\"

CXXFLAGS := -std=$(STD) $(WARN)
CFLAGS   ?= -O3 -DNDEBUG
LDFLAGS  :=
LDLIBS   := -lfftw3 -lvoro++ -L./src/voro++-0.4.6/src

ifeq ($(OMP),1)
  CXXFLAGS += -fopenmp
  LDLIBS   += -fopenmp
endif

# Build profile
ifeq ($(MODE),debug)
  CXXFLAGS += -O0 -g3
else ifeq ($(MODE),release)
  CXXFLAGS += -O3 -DNDEBUG
else
  $(error Unknown MODE '$(MODE)'; use MODE=debug or MODE=release)
endif

# Optional: ccache if installed
ifneq ($(shell command -v ccache 2>/dev/null),)
  CXX := ccache $(CXX)
endif

# --- xdrfile (GROMACS XTC/TRR) integration ----------------------------------
# The bundled C library lives at ./xdrfile-1.1b
XDR_DIR      ?= xdrfile-1.1b
XDR_INC       = $(XDR_DIR)/include
XDR_SRC_DIR   = $(XDR_DIR)/src

# Make sure C and C++ compiles can find headers
CPPFLAGS += -I./$(XDR_INC) -I./$(XDR_DIR)
XDR_CFLAGS := -I./$(XDR_INC) -I./$(XDR_DIR) $(CFLAGS)

# C sources to build from xdrfile
XDR_CSRCS := $(XDR_SRC_DIR)/xdrfile.c \
             $(XDR_SRC_DIR)/xdrfile_xtc.c \
             $(XDR_SRC_DIR)/xdrfile_trr.c
XDR_OBJS  := $(patsubst $(XDR_SRC_DIR)/%.c,$(BUILD_DIR)/xdr/%.o,$(XDR_CSRCS))

# Some distros require -lm for math symbols used by xdrfile
LDLIBS += -lm

# Ensure dir exists
$(shell mkdir -p $(BUILD_DIR)/xdr)

# --- Source discovery (flat, only src/*.cpp) ---------------------------------
SRCS := $(wildcard $(SRC_DIR)/*.cpp)

# Exclude specific sources
EXCLUDE_SRCS := src/time_correlation_function.cpp \
                src/spherical_wave_vectors.cpp \
								src/spacial_decomposition.cpp
SRCS := $(filter-out $(EXCLUDE_SRCS),$(SRCS))

# Map to objects in build/
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

# Ensure build dirs exist
$(shell mkdir -p $(BUILD_DIR))

# --- Targets -----------------------------------------------------------------
.PHONY: all clean distclean rebuild format lint help
all: $(APP)

.PHONY: voro
voro:
	$(MAKE) -C src/voro++-0.4.6

# Final link: include xdrfile objects as well
$(APP): $(OBJS) $(XDR_OBJS) $(VORO_OBJS) | voro
	@echo "  LINK    $@"
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Compile C++: strictly src/<file>.cpp → build/<file>.o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "  CXX     $<"
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Compile C from xdrfile: xdrfile-1.1b/src/*.c → build/xdr/*.o
$(BUILD_DIR)/xdr/%.o: $(XDR_SRC_DIR)/%.c
	@echo "  CC      $<"
	$(CC) $(XDR_CFLAGS) -c $< -o $@

## Compile voro++: src/voro++-0.4.6/src/*.cc → build/voro/*.o (isolated includes)
#$(BUILD_DIR)/voro/%.o: $(VORO_DIR)/%.cc
#	@echo "  CXX     $< (voro++)"
#	$(CXX) $(VORO_CPPFLAGS) $(CXXFLAGS) -c $< -o $@
#
# Include auto-generated header dependencies (for C++ objects)
-include $(DEPS)

clean:
	@echo "  CLEAN   objects"
	@rm -rf $(BUILD_DIR)/*

distclean: clean
	@echo "  CLEAN   binary"
	@rm -f $(APP)

rebuild: distclean all

# Code quality helpers (optional)
format:
	@command -v clang-format >/dev/null || { echo "clang-format not found"; exit 1; }
	@echo "  FORMAT  src/**/*.h,*.hpp,*.cpp"
	@find src -type f \( -name '*.hpp' -o -name '*.h' -o -name '*.cpp' \) -print0 | xargs -0 clang-format -i

lint:
	@command -v clang-tidy >/dev/null || { echo "clang-tidy not found"; exit 1; }
	@echo "  LINT    (clang-tidy)"
	@clang-tidy -p . $(SRCS)

help:
	@echo "Targets:"
	@echo "  all (default)  Build $(APP) [MODE=$(MODE) OMP=$(OMP)]"
	@echo "  clean          Remove objects only"
	@echo "  distclean      Remove objects and binary"
	@echo "  rebuild        Full clean + build"
	@echo "  format         Run clang-format (if installed)"
	@echo "  lint           Run clang-tidy (if installed)"
	@echo
	@echo "Examples:"
	@echo "  make MODE=debug"
	@echo "  make MODE=release OMP=0"
