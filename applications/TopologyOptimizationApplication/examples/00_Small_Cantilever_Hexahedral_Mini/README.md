# 00 - Small Cantilever Hexahedral - Mini
[description]
modeled by Philipp Hofer.

## Model parameters

### Unit system used
[unit system]

### Geometry
[geometrical description]


### Material
Steel

Young's modulus: 

Poison's ratio:

Density:

### Solver settings
Solution type: static analysis 

### Model details
Number of nodes: 6408

Number of elements: 4800

Element type: Small displacement volume - hexahedron

### Figures
[Problem]

[Optimal topology]

[Smoothed]


### Notes
geändert:

1.	runn_TopOpt.py:
		import hinzugefügt
		solver_module anstatt import_module

2. 	ProjectParameter:
		StaticSIMP als solvertype auf Static umgeändert

Versuche es mit der Vorlage aus GitHub das run_TopOpt zusammenzustellen. Bei Versuchen schon weiter gekommen als mit momentaner Variante!

Aufpassen mit Solver. Muss dieser selbst erstellt werden? Oder wird ein bereits vorhanderner Solver verwendet? Ist es das was im Programm steht und warum es nicht weiter geht? Den neuen Solver in python_scripts nochmal anschauen und vergleichen mit den Solver für SructuralMechanicsApplication


