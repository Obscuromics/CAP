install:
	@echo "Setting up Conda environment..."
	bash setup_conda.sh
	@echo "Compiling BCT C++ alternative..."
	# We try to use the compiler from the new environment if possible, or system compiler
	# This assumes 'conda run' works
	conda run -n cap-pipeline make -C bin/src/BCT

compile:
	make -C bin/src/BCT

clean:
	@echo "Removing Conda environment..."
	conda env remove -n cap-pipeline
	rm -f bin/ctw-calc
