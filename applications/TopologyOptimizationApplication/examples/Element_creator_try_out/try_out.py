import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as ksm
import KratosMultiphysics.TopologyOptimizationApplication as kto



current_model = km.Model()
mp = current_model.CreateModelPart("solid_part")

mp.CreateNewNode(1, 0.00000,  1.00000,  1.00000)
mp.CreateNewNode(2, 0.16500,  0.74500,  0.70200)
mp.CreateNewNode(3, 0.27300,  0.75000,  0.23000)
mp.CreateNewNode(4, 0.78800,  0.69300,  0.64400)
mp.CreateNewNode(5, 0.32000,  0.18600,  0.64300)
mp.CreateNewNode(6, 0.00000,  1.00000,  0.00000)
mp.CreateNewNode(7, 0.00000,  0.00000,  1.00000)
mp.CreateNewNode(8, 1.00000,  1.00000,  1.00000)
mp.CreateNewNode(9, 0.67700,  0.30500,  0.68300)
mp.CreateNewNode(10, 0.24900,  0.34200,  0.19200)
mp.CreateNewNode(11, 0.85000,  0.64900,  0.26300)
mp.CreateNewNode(12, 0.82600,  0.28800,  0.28800)
mp.CreateNewNode(13, 0.00000,  0.00000,  0.00000)
mp.CreateNewNode(14, 1.00000,  1.00000,  0.00000)
mp.CreateNewNode(15, 1.00000,  0.00000,  1.00000)
mp.CreateNewNode(16, 1.00000,  0.00000,  0.00000)

#create Element
mp.CreateNewElement("SmallDisplacementSIMPElement3D8N", 1,[10,5,2,3,13,7,1,6], mp.GetProperties()[1])
mp.CreateNewElement("SmallDisplacementSIMPElement3D8N", 2,[12,9,5,10,16,15,7,13], mp.GetProperties()[1])
mp.CreateNewElement("SmallDisplacementSIMPElement3D8N", 3,[12,11,3,10,9,4,2,5], mp.GetProperties()[1])
mp.CreateNewElement("SmallDisplacementSIMPElement3D8N", 4,[9,4,2,5,15,8,1,7], mp.GetProperties()[1])
mp.CreateNewElement("SmallDisplacementSIMPElement3D8N", 5,[4,11,3,2,8,14,6,1], mp.GetProperties()[1])
mp.CreateNewElement("SmallDisplacementSIMPElement3D8N", 6,[11,4,9,12,14,8,15,16], mp.GetProperties()[1])
mp.CreateNewElement("SmallDisplacementSIMPElement3D8N", 7,[11,12,10,3,14,16,13,6], mp.GetProperties()[1])
print("Done")
