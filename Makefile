.PHONY: clean

objects :=\
	plots/troubleshootAPA_before.pdf\
	plots/troubleshootAPA_after.pdf

all: $(objects)

clean:
	rm -rf $(objects)

plots/troubleshootAPA_before.pdf\
plots/troubleshootAPA_after.pdf:\
	scripts/utils/findBadLoops.R\
	scripts/utils/layoutFunctions.R\
	scripts/utils/plotApas.R\
	data/raw/hic/hicFiles/YAPP/YAPP_HEK_control_inter_30.hic\
	data/raw/hic/hicFiles/YAPP/YAPP_HEK_sorbitol_inter_30.hic\
	data/raw/hic/loops/YAPP_hic_diff_loopCounts.rds\
	scripts/troubleshootAPA.R
		mkdir -p plots
		Rscript scripts/troubleshootAPA.R
		