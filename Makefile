install:
	@echo "Setting up Conda environment..."
	bash setup_conda.sh

compile:
	make -C bin/src/BCT

clean:
	@echo "Removing Conda environment..."
	conda env remove -n cap-pipeline
	rm -f bin/ctw-calc
