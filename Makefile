install:
	@echo "Setting up Conda environment..."
	bash setup_conda.sh
	@echo "Compiling BCT C++ alternative..."
	# We try to use the compiler from the new environment if possible, or system compiler
	# The '-' prefix makes this step non-fatal (ignores errors)
	-conda run -n cap-pipeline make -C bin/src/BCT || echo "⚠️  C++ compilation failed. The pipeline will rely on the 'BCT' R package."
	@echo "Setting permissions..."
	chmod +x modules/TRASH_2/src/TRASH.R

compile:
	make -C bin/src/BCT

clean:
	@echo "Removing Conda environment..."
	conda env remove -n cap-pipeline
	rm -f bin/ctw-calc
