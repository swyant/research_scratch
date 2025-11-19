import dakota.interfacing as di
import subprocess 
print("HEYYY")
params, results = di.read_parameters_file()

di.dprepro(template="plate_with_hole_input.template", 
           parameters=params, 
           output="plate_with_hole_input")

subprocess.run(["uv", "run", "python", "plate_with_hole.py"])

## post proces 
with open("plate_with_hole_results.txt", "r") as f: 
    max_vm_stress = float(f.readline().strip())

results['max_von_Mises_stres'].function = max_vm_stress
results.write()
