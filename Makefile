CASES := $(wildcard cases/*)
resfiles := $(foreach dir,$(CASES),$(dir)/res.pdf)
resfiles2 := $(foreach dir,$(CASES),$(dir)/res.tab)
resfiles3 := $(foreach dir,$(CASES),$(dir)/sdrep.tab)
resfiles4 := $(foreach dir,$(CASES),$(dir)/jit.tab)

datfiles := $(foreach dir,$(CASES),$(dir)/Y.tab)

.PHONY: sim clean run printcases allrestab summary compile

cases/sim1/Y.tab: src/sim1.R
	echo 'source("src/sim1.R")' | R --slave

sim: cases/sim1/Y.tab

clean:
	rm -f $(resfiles) $(resfiles2) $(resfiles3)  $(resfiles4) src/*.o src/*.so sumtab.txt meantab.txt

src/gmrf1.so: src/gmrf1.cpp
	echo 'TMB:::compile("src/gmrf1.cpp")' | R --slave

compile: src/gmrf1.so

res.pdf: Y.tab ../../src/gmrf1.so ../../src/run.R 
	echo 'source("../../src/run.R")' | R --slave

$(resfiles): src/gmrf1.so src/run.R $(datfiles)
	@$(MAKE) -C $(@D) -i -s -f ../../Makefile res.pdf

run: $(resfiles)

printcases:
	$(info $$CASES are [$(CASES)])
	$(info $$resfiles are [$(resfiles)])
	$(info $$datfiles are [$(datfiles)])

allrestab:
	cat $(resfiles2) > allrestab.txt

sumtab.txt: $(resfiles)
	cd src; echo 'source("summarize.R")' | R --slave

summary: sumtab.txt 
