TOOLS_BINARIES = tools/rmapAnaly/rmap-analy \
	tools/CompToPCAlign/comp2pcalign \
	tools/scripts/merge-fastq.py \
	tools/scripts/count_errors.py \
	tools/scripts/ecutils.py \
	tools/scripts/hitec-analy.pl \
	tools/scripts/compute-stats.py \
	tools/scripts/fastq-converter-v2.0.pl \
	tools/scripts/coral-analy.pl \
	tools/scripts/quake-analy.py \
	tools/scripts/alignment.py \
	tools/scripts/sam-analysis.py

TOOLS_ZIP_DATA = tools/rmapAnaly/Errlocator.h \
	tools/rmapAnaly/rmapAnaly.h \
	tools/rmapAnaly/main.cpp \
	tools/rmapAnaly/fasta_file.hpp \
	tools/rmapAnaly/Errlocator.cpp \
	tools/rmapAnaly/makefile \
	tools/scripts \
	tools/scripts/merge-fastq.py \
	tools/scripts/count_errors.py \
	tools/scripts/ecutils.py \
	tools/scripts/hitec-analy.pl \
	tools/scripts/compute-stats.py \
	tools/scripts/fastq-converter-v2.0.pl \
	tools/scripts/coral-analy.pl \
	tools/scripts/quake-analy.py \
	tools/scripts/alignment.py \
	tools/scripts/sam-analysis.py \
	tools/CompToPCAlign \
	tools/CompToPCAlign/Record.cpp \
	tools/CompToPCAlign/util.h \
	tools/CompToPCAlign/main.cpp \
	tools/CompToPCAlign/compSR.h \
	tools/CompToPCAlign/Record.h \
	tools/CompToPCAlign/makefile

OUTPUT_ORG = doc/ecet.html

help:
	@echo 'Options available : '
	@echo 'make zip           make a zip file of all the source code'
	@echo 'make ecet          build all the executables and copy executables/scripts to bin directory'
	@echo 'make doc           generate doc : needs emacs >= 25.1 org mode 7'
	@echo 'make clean         remove all built files, documentation and zip'
	@echo ''

doc: $(OUTPUT_ORG)

$(OUTPUT_ORG): %.html : %.org
	emacs --batch --visit=$< --execute='(org-mode)' --execute='(org-html-export-to-html)'


ecet:
	make -C tools/rmapAnaly
	make -C tools/CompToPCAlign
	mkdir -p bin
	cp $(TOOLS_BINARIES) bin/

clean:
	make -C tools/rmapAnaly clean
	make -C tools/CompToPCAlign clean
	-rm -rf bin
	-rm tools.zip

zip:
	-rm tools.zip
	zip tools.zip $(TOOLS_ZIP_DATA)
