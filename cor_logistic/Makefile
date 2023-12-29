## all           : clean saved results and get results to populate report
all: clean cor_logistic report

report:
	bash _build_chapter.sh


## clean         : delete saved bootstrap results
clean:
	rm -rf output/$(TRIAL)/*

## cor_logistic     : create tables and figures
cor_logistic:
ifeq ($(TRIAL),$(filter $(TRIAL), id27hpv))
	Rscript code/cor_logistic.R M18
	Rscript code/cor_logistic.R M18sus
	Rscript code/cor_logistic.R M18nAb
	Rscript code/cor_logistic.R M18susnAb
endif


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_logistic
