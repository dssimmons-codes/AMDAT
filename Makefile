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

# Wavevector bootstrap
WV_STAMP := $(SRC_DIR)/qvectors/.ready

# Add version.h for automated version info at build time
VERSION_H := src/generated/version.h

# Try to get a friendly describe string; ignore errors if not a repo
GIT_DESCRIBE := $(shell git describe --tags --dirty --always 2>/dev/null)
GIT_COMMIT   := $(shell git rev-parse --short HEAD 2>/dev/null)
GIT_BRANCH   := $(shell git rev-parse --abbrev-ref HEAD 2>/dev/null)

# Extract SemVer from tag if present (v1.2.3 → 1.2.3)
# Fallback: read VERSION, else default
ifeq ($(strip $(GIT_DESCRIBE)),)
  AMDAT_SEMVER := $(strip $(shell [ -f VERSION ] && cat VERSION))
  ifeq ($(strip $(AMDAT_SEMVER)),)
    AMDAT_SEMVER := 0.0.0+archive
  endif
else
  # pull the first thing that looks like vX.Y.Z or X.Y.Z
  AMDAT_SEMVER := $(shell echo "$(GIT_DESCRIBE)" | sed -E 's/.*(v?[0-9]+\.[0-9]+\.[0-9]+).*/\1/' | sed 's/^v//')
endif

# Allow packagers to override from environment: AMDAT_SEMVER=1.0.0+pkg make
ifneq ($(origin AMDAT_SEMVER),file)
  AMDAT_SEMVER := $(AMDAT_SEMVER)
endif

BUILD_DATE   := $(shell date -u +%Y-%m-%dT%H:%M:%SZ)

# --- Flags -------------------------------------------------------------------
CPPFLAGS := -MMD -MP -I./src \
						-I./src/generated \
						-I./third_party/voro++-0.4.6/src \
						-L./third_party/voro++-0.4.6/src \
						-I./third_party/xdrfile-1.1b/src \
            -DWV3D=\"$(WV3D)\" \
            -DWV2D=\"$(WV2D)\" \
            -DWV1D=\"$(WV1D)\"

CXXFLAGS := -std=$(STD) $(WARN)
CFLAGS   ?= -O3 -DNDEBUG
LDFLAGS  :=
LDLIBS   := -lvoro++ -L./third_party/voro++-0.4.6/src
FFTW_ROOT ?=

# --- try pkg-config first ---
PKGCONF  ?= pkg-config
ifeq ($(shell which $(PKGCONF) >/dev/null 2>&1 && echo yes),yes)
  ifeq ($(shell $(PKGCONF) --exists fftw3 && echo yes),yes)
    CPPFLAGS += $(shell $(PKGCONF) --cflags fftw3)
    LDLIBS   += $(shell $(PKGCONF) --libs   fftw3)
  endif
endif

# --- fallback: user-specified prefix ---
ifneq ($(FFTW_ROOT),)
  CPPFLAGS += -I$(FFTW_ROOT)/include
  LDFLAGS  += -L$(FFTW_ROOT)/lib
  LDLIBS   += -lfftw3
endif

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
# The bundled C library lives at ./third_party/xdrfile-1.1b
XDR_DIR      ?= third_party/xdrfile-1.1b
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

$(VERSION_H):
	@mkdir -p $(dir $@)
	@{ \
    echo '#pragma once'; \
    echo 'namespace amdat::build {'; \
    echo 'inline constexpr const char* NAME="AMDAT";'; \
    echo 'inline constexpr const char* SEMVER="$(AMDAT_SEMVER)";'; \
    echo 'inline constexpr const char* GIT_DESCRIBE="$(GIT_DESCRIBE)";'; \
    echo 'inline constexpr const char* GIT_COMMIT="$(GIT_COMMIT)";'; \
    echo 'inline constexpr const char* GIT_BRANCH="$(GIT_BRANCH)";'; \
    echo 'inline constexpr const char* BUILD_DATE="$(BUILD_DATE)";'; \
    echo 'inline constexpr const char* COMPILER='; \
    echo '#if defined(__clang__)'; \
    echo '  "clang " __clang_version__;'; \
    echo '#elif defined(__GNUC__)'; \
    echo '  "gcc " __VERSION__;'; \
    echo '#elif defined(_MSC_VER)'; \
    echo '  "MSVC";'; \
    echo '#else'; \
    echo '  "unknown";'; \
    echo '#endif'; \
    echo 'inline constexpr bool OPENMP='; \
    echo '#ifdef _OPENMP'; \
    echo 'true;'; \
    echo '#else'; \
    echo 'false;'; \
    echo '#endif'; \
    echo '}'; \
  } > $@


.PHONY: voro
voro:
	$(MAKE) -C third_party/voro++-0.4.6 \
					CC="$(CC)" CXX="$(CXX)"

.PHONY: clean_voro
clean_voro:
	$(MAKE) -C third_party/voro++-0.4.6 clean

.PHONY: clean_qvectors
clean_qvectors:
	@$(RM) -rf $(SRC_DIR)/qvectors/.ready
	@$(MAKE) -C $(SRC_DIR)/qvectors distclean

.PHONY: qvectors
qvectors: $(WV_STAMP)

$(WV_STAMP):
	@$(MAKE) -C $(SRC_DIR)/qvectors

# Final link: include xdrfile objects as well
$(APP): $(OBJS) $(XDR_OBJS) $(VORO_OBJS) | qvectors voro
	@echo "  LINK    $@"
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Compile C++: strictly src/<file>.cpp → build/<file>.o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(VERSION_H)
	@echo "  CXX     $<"
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Compile C from xdrfile: third_party/xdrfile-1.1b/src/*.c → build/xdr/*.o
$(BUILD_DIR)/xdr/%.o: $(XDR_SRC_DIR)/%.c
	@echo "  CC      $<"
	$(CC) $(XDR_CFLAGS) -c $< -o $@

## Compile voro++: third_party/voro++-0.4.6/src/*.cc → build/voro/*.o (isolated includes)
#$(BUILD_DIR)/voro/%.o: $(VORO_DIR)/%.cc
#	@echo "  CXX     $< (voro++)"
#	$(CXX) $(VORO_CPPFLAGS) $(CXXFLAGS) -c $< -o $@
#
# Include auto-generated header dependencies (for C++ objects)
-include $(DEPS)

clean: clean_voro
	@echo "  CLEAN   objects"
	@rm -rf $(BUILD_DIR)/*

distclean: clean clean_qvectors
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

# --- Conda helpers (do not require 'conda activate') -------------------------
CONDA       ?= conda
ENV_NAME    ?= amdat
# Prefer an exact lock if present, else fall back to environment.yml
ENV_SPEC    ?= $(if $(wildcard conda-linux-64.txt),conda-linux-64.txt,environment.yml)
ENV_STAMP   := .env/conda.ready

# Create/update the Conda environment
.PHONY: conda-setup
conda-setup: $(ENV_STAMP)

$(ENV_STAMP):
	@echo ">> Creating conda env '$(ENV_NAME)' from $(ENV_SPEC)"
	@mkdir -p .env
	@if [ "$(ENV_SPEC)" = "conda-linux-64.txt" ]; then \
	  $(CONDA) create -y -n $(ENV_NAME) --file $(ENV_SPEC); \
	else \
	  $(CONDA) env create -n $(ENV_NAME) -f $(ENV_SPEC) || $(CONDA) env update -n $(ENV_NAME) -f $(ENV_SPEC); \
	fi
	@touch $@

# Build *inside* the Conda env, ensuring:
#  - we use conda-forge compilers if present
#  - pkg-config sees FFTW from the env
#  - parallel build matches CPU count on the node
.PHONY: conda
conda: $(ENV_STAMP)
	@$(CONDA) run -n $(ENV_NAME) bash -lc '\
	  set -euo pipefail; \
	  : "$${CONDA_PREFIX:?conda run did not set CONDA_PREFIX}"; \
	  export PKG_CONFIG_PATH="$$CONDA_PREFIX/lib/pkgconfig:$${PKG_CONFIG_PATH:-}"; \
	  # Prefer conda-forge triplet compilers if available; otherwise keep Makefile defaults \
	  if command -v x86_64-conda-linux-gnu-g++ >/dev/null 2>&1; then \
	    MAKE_CXX="CXX=x86_64-conda-linux-gnu-g++"; \
	    MAKE_CC="CC=x86_64-conda-linux-gnu-gcc"; \
	  else \
	    MAKE_CXX=""; MAKE_CC=""; \
	  fi; \
	  echo ">> Using $$MAKE_CC $$MAKE_CXX"; \
	  export CPATH="$$CONDA_PREFIX/include:$${CPATH:-}"; \
	  export LIBRARY_PATH="$$CONDA_PREFIX/lib:$${LIBRARY_PATH:-}"; \
		echo ">> pkg-config: $(which pkg-config || true)"; \
		echo ">> fftw cflags: $$(pkg-config --cflags fftw3 || true)"; \
		echo ">> fftw libs:   $$(pkg-config --libs   fftw3 || true)"; \
	  $(MAKE) $$MAKE_CC $$MAKE_CXX MODE=$(MODE) OMP=$(OMP) FFTW_ROOT="$$CONDA_PREFIX" \
	'

# Optional: interactive shell in the env (handy for debugging)
.PHONY: conda-shell
conda-shell: $(ENV_STAMP)
	@$(CONDA) run -n $(ENV_NAME) bash

# Clean just your build outputs (keeps the env)
.PHONY: conda-clean
conda-clean:
	@$(RM) -rf $(BUILD_DIR)/*

# Nukes the env stamp (does NOT delete the conda env itself)
.PHONY: conda-unstamp
conda-unstamp:
	@$(RM) -f $(ENV_STAMP)
