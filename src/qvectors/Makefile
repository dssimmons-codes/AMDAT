# Unpacks qvectors.zip into this directory, once, producing a stamp file
# in src/qvectors/ and qvectors1d, qvectors2d, and qvectors3d inside it.
QZIP       := qvectors.zip      # committed to git
REQ_DIRS   := qvectors1d qvectors2d qvectors3d
REQ_STAMPS := $(addsuffix /.dir,$(REQ_DIRS))
QSTMP      := .ready            # stamp to mark successful unpack

QZIP := $(strip $(QZIP))
QSTMP := $(strip $(QSTMP))

.PHONY: all clean distclean
all: $(QSTMP)

# Master stamp depends on per-dir stamps
$(QSTMP): $(REQ_STAMPS)
	@touch "$(QSTMP)"

# We unpack once, then (re)create all per-dir stamps.
$(REQ_STAMPS): $(QZIP)
	@echo ">> Unzipping $(QZIP) into $(CURDIR)"
	@unzip -q -o "$(QZIP)"
	@touch "$(QSTMP)"
	@for d in $(REQ_DIRS); do \
		test -d "$$d" || { echo "!! Missing $$d after unzip"; exit 1; }; \
	done
	@# refresh stamps for all required dirs
	@for d in $(REQ_DIRS); do \
		touch "$$d/.dir"; \
	done

clean:
	@rm -rf "$(QSTMP)" $(REQ_DIRS)

distclean: clean
	@echo ">> Removing unpacked qvectors"
	@find . -mindepth 1 -maxdepth 1 \
	  ! -name "$(notdir $(QZIP))" \
	  ! -name 'Makefile' \
	  ! -name '.gitignore' \
		! -name '*.cpp' \
		! -name '*.m' \
	  -exec rm -rf {} +
