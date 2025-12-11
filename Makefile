install:
	@echo "Setting up Conda environment..."
	bash setup_conda.sh

clean:
	@echo "Removing Conda environment..."
	conda env remove -n cap-pipeline
