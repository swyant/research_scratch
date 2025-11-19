import felupe as fem
import numpy as np
## geometry ##
h = 1
L = 2
r = 0.3
mesh_file = "partial_plate_with_hole.vtk" 
####

## constitutive params ## 
# using a separate file so that Dakota can do the param interpolation
# without running into {...} conflicts with the string.format below
with open("plate_with_hole_input", "r") as f:
  E = float(f.readline().strip())
  nu = float(f.readline().strip())
#### 

mesh = fem.mesh.read("partial_plate_with_hole.vtk",dim=2)[0]
region = fem.RegionQuad(mesh)
displacement = fem.Field(region, dim=2)
field = fem.FieldContainer([displacement])

boundaries = fem.dof.symmetry(displacement)

umat = fem.LinearElasticPlaneStress(E=210000, nu=0.3)
solid = fem.SolidBody(umat, field)

sig_inf=-100
region_boundary = fem.RegionQuadBoundary(mesh, mask=mesh.points[:, 0] == L)
field_boundary = fem.FieldContainer([fem.Field(region_boundary, dim=2)])

load = fem.SolidBodyPressure(field_boundary, pressure=sig_inf)

step = fem.Step(items=[solid, load], boundaries=boundaries)
job = fem.Job(steps=[step]).evaluate()

max_vm_stress = float(np.max(fem.tools.project(solid.evaluate.stress(), region)))

with open("plate_with_hole_results.txt", "w") as f:
    f.write("{:.18f}\n".format(max_vm_stress))
