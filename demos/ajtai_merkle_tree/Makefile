# Convenience targets for the Ajtai-Merkle / LaBRADOR demo.
# Run `conda activate sage` first (provides python3 + the lazer module + sage).
.PHONY: all l2 soundness estimate clean

all: l2 soundness

l2:
	python3 example_l2norm.py

soundness:
	python3 test_soundness.py

estimate:
	sage estimate_sis_params.sage

clean:
	rm -rf __pycache__ estimate_sis_params.sage.py
