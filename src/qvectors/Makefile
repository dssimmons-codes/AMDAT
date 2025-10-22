# Unpacks qvectors.zip into this directory, once, producing a stamp file.

QZIP  := qvectors.zip      # committed to git
QSTMP := .ready            # stamp to mark successful unpack

QZIP := $(strip $(QZIP))
QSTMP := $(strip $(QSTMP))

.PHONY: all clean distclean
all: $(QSTMP)

$(QSTMP): $(QZIP)
	@echo ">> Unzipping $(QZIP) into $(CURDIR)"
	@unzip -q -o "$(QZIP)"
	@touch "$(QSTMP)"

# Remove generated payload but keep the archive and this Makefile
clean:
	@rm -f "$(QSTMP)"

# If you want a hard reset of the directory (except the zip/Makefile)
distclean: clean
	@echo ">> Removing unpacked qvectors"
	@find . -mindepth 1 -maxdepth 1 \
	  ! -name "$(notdir $(QZIP))" \
	  ! -name 'Makefile' \
	  ! -name '.gitignore' \
		! -name '*.cpp' \
		! -name '*.m' \
	  -exec rm -rf {} +
