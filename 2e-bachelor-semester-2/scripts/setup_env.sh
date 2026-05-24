#!/usr/bin/env bash
set -euo pipefail

# Create and activate a Python virtualenv named .venv, then install Python deps
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install -r requirements.txt

# Setup Julia project for Gmsh examples (requires Julia installed)
if command -v julia >/dev/null 2>&1; then
  julia -e 'import Pkg; Pkg.activate(joinpath(pwd(),"julia")); Pkg.add("Gmsh"); Pkg.instantiate()'
else
  echo "Julia not found. Install Julia and then run:"
  echo "  julia -e 'import Pkg; Pkg.activate(joinpath(pwd(),"\"julia\"")); Pkg.add(\"Gmsh\"); Pkg.instantiate()'"
fi

echo "Setup complete. Activate the venv with: source .venv/bin/activate"
