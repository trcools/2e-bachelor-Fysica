# Compatibility wrapper for Julia examples using Gmsh
try
  import Gmsh
  const gmsh = Gmsh.gmsh
catch e
  @error "Gmsh package not found; run 'import Pkg; Pkg.add(\"Gmsh\")' in Julia" error=e
  rethrow()
end
