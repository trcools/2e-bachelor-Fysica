Setup for gmsh / examples
=========================

Python (recommended)
- Create and activate the virtual environment and install Python deps:

```bash
./setup_env.sh
source .venv/bin/activate
```

This will install `numpy` and the `gmsh` Python wheel into `.venv`.

Run Python examples using the venv python, for example:

```bash
${workspaceFolder}/Sterrenstelsels/.venv/bin/python ${workspaceFolder}/.venv/share/doc/gmsh/examples/api/adapt_mesh.py -nopopup
```

VS Code
- The workspace setting `.vscode/settings.json` points to the recommended interpreter:
[.vscode/settings.json](.vscode/settings.json#L1)

Julia examples
- Prepare the Julia project and add the `Gmsh` package:

```bash
julia -e 'import Pkg; Pkg.activate(joinpath(pwd(),"julia")); Pkg.add("Gmsh"); Pkg.instantiate()'
```

- Use the compatibility wrapper `julia/gmsh_compat.jl` from Julia examples to access the `gmsh` API:

```julia
include("julia/gmsh_compat.jl")
# then normal example code can use `gmsh` as in the upstream examples
```
