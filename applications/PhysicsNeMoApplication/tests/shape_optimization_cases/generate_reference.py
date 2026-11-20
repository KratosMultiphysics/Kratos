"""Reference vertex-morphing fields from ShapeOptimizationApplication.

Run in a WHEEL-ONLY environment (KratosMultiphysics + KratosShapeOptimization
from PyPI). Never mixed with a locally compiled core: the wheels are GCC-built
and the local core is Clang-built, and pybind11 keys its type registry on
compiler identity.
"""
import numpy, KratosMultiphysics as Kratos
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory

DIV, RADIUS = 6, 0.4
FILTERS = ("linear", "gaussian")

def BuildSurface(model, name):
    mp = model.CreateModelPart(name)
    mp.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
    for v in (KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE, Kratos.NORMAL,
              KSO.NORMALIZED_SURFACE_NORMAL, KSO.SHAPE_CHANGE):
        mp.AddNodalSolutionStepVariable(v)
    props = mp.CreateNewProperties(1)
    nid = lambda i, j: i * (DIV + 1) + j + 1
    for i in range(DIV + 1):
        for j in range(DIV + 1):
            mp.CreateNewNode(nid(i, j), i / DIV, j / DIV, 0.0)
    e = 0
    for i in range(DIV):
        for j in range(DIV):
            e += 1; mp.CreateNewElement("Element3D3N", e, [nid(i,j), nid(i+1,j), nid(i+1,j+1)], props)
            e += 1; mp.CreateNewElement("Element3D3N", e, [nid(i,j), nid(i+1,j+1), nid(i,j+1)], props)
    return mp

out = {}
model = Kratos.Model()
mp = BuildSurface(model, "Design")
coords = numpy.array([[n.X, n.Y, n.Z] for n in mp.Nodes])
ids = numpy.array([n.Id for n in mp.Nodes], dtype=numpy.int64)
out["coordinates"] = coords
out["node_ids"] = ids
out["radius"] = numpy.array([RADIUS])

# two control cases: a single impulse, and a uniform field (translation test)
centre_id = int(ids[len(ids)//2])
cases = {"impulse": {centre_id: [0.0, 0.0, 1.0]},
         "uniform": {int(i): [0.0, 0.0, 1.0] for i in ids}}

for filter_name in FILTERS:
    for case_name, control in cases.items():
        model2 = Kratos.Model()
        part = BuildSurface(model2, "Design")
        settings = Kratos.Parameters("""{
            "filter_function_type"       : "linear",
            "filter_radius"              : 0.4,
            "max_nodes_in_filter_radius" : 10000,
            "matrix_free_filtering"       : false,
            "consistent_mapping"          : false,
            "improved_integration"        : false
        }""")
        settings["filter_function_type"].SetString(filter_name)
        settings["filter_radius"].SetDouble(RADIUS)
        mapper = mapper_factory.CreateMapper(part, part, settings)
        mapper.Initialize()
        for node in part.Nodes:
            node.SetSolutionStepValue(KSO.CONTROL_POINT_UPDATE,
                                      control.get(node.Id, [0.0, 0.0, 0.0]))
        mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
        field = numpy.array([list(n.GetSolutionStepValue(KSO.SHAPE_UPDATE))
                             for n in part.Nodes])
        out[f"{filter_name}_{case_name}"] = field
        print(f"R: {filter_name}/{case_name}: |field|max={numpy.abs(field).max():.6f} "
              f"nonzero={int((numpy.abs(field).sum(axis=1) > 1e-12).sum())}/{len(field)}")

out["centre_node_id"] = numpy.array([centre_id], dtype=numpy.int64)
numpy.savez_compressed("/tmp/ksovenv/vertex_morphing_reference.npz", **out)
print("R: saved keys:", sorted(out))
