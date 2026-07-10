all: preprocess rda

preprocess:
	Rscript run_pipeline.R --step preprocess

rda:
	Rscript run_pipeline.R --step rda

smoke_test:
	Rscript data/sample/generate_sample.R
	Rscript run_pipeline.R --step preprocess
	Rscript run_pipeline.R --step rda

clean:
	rm -rf data/processed/* results/*

.PHONY: all preprocess rda smoke_test clean
