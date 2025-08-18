## all           : clean saved results and get results to populate report
all: clean cor_poisson report

report:
	bash _build_chapter.sh


## clean         : delete saved bootstrap results
clean:
	rm -rf output/$(TRIAL)/*

## cor_poisson     : create tables and figures
cor_poisson:
ifeq ($(TRIAL),$(filter $(TRIAL), iliad_ib202p))
	Rscript code/cor_poisson_iliad_ib202p.R C0iliad_ib202p
endif


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_poisson
