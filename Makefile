CASES := $(wildcard testmore/*/)

.PHONY: sim clean runsim

cases/sim1/Y.tab: src/sim1.R
	echo 'source("src/sim1.R")' | R --slave

sim: cases/sim1/Y.tab

clean:
	rm -f cases/sim1/res.pdf src/*.o src/*.so 

cases/sim1/res.pdf: src/run.R cases/sim1/Y.tab
	cd cases/sim1; echo 'source("../../src/run.R")' | R --slave

runsim: cases/sim1/res.pdf

