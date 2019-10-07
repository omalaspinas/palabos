SUBDIRS := $(wildcard */.)
SUBDIRS := $(filter-out include/., $(SUBDIRS))
SUBDIRSCLEAN=$(addsuffix \ clean ,$(SUBDIRS))

.PHONY : all $(SUBDIRS) clean $(SUBDIRSCLEAN)

all : $(SUBDIRS)
clean: $(SUBDIRSCLEAN)

$(SUBDIRS) :
	$(MAKE) -C $@

$(SUBDIRSCLEAN):
	$(MAKE) -C $@