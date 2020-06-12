CASES := $(wildcard cases/*)
resfiles := $(foreach dir,$(CASES),$(dir)/res.pdf)
datfiles := $(foreach dir,$(CASES),$(dir)/Y.tab)

.PHONY: sim clean run printcases

cases/sim1/Y.tab: src/sim1.R
	echo 'source("src/sim1.R")' | R --slave

sim: cases/sim1/Y.tab

clean:
	rm -f $(resfiles) src/*.o src/*.so

res.pdf: Y.tab ../../src/gmrf1.cpp ../../src/run.R
	echo 'source("../../src/run.R")' | R --slave

$(resfiles): src/gmrf1.cpp src/run.R $(datfiles)
	@$(MAKE) -C $(@D) -i -s -f ../../Makefile res.pdf

run: $(resfiles)

printcases:
	$(info $$CASES are [$(CASES)])
	$(info $$resfiles are [$(resfiles)])
	$(info $$datfiles are [$(datfiles)])
