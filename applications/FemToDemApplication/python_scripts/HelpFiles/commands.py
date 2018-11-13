
for node in self.DEM_Solution.spheres_model_part.Nodes:
	#print(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
	print(node.GetSolutionStepValue(KratosMultiphysics.TOTAL_FORCES_X))
	#print(node.Id)
	#node.SetSolutionStepValue(KratosMultiphysics.RADIUS, 0.1)


for node in self.FEM_Solution.main_model_part.Nodes:
print(node.GetSolutionStepValue(KratosMultiphysics.TOTAL_FORCES_X))


#print("esferas" , self.SpheresModelPart)
#Wait()
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

# test with two balls to test the DEM Application
Coordinates1 = KratosMultiphysics.Array3()
Coordinates1[0] = 0.0
Coordinates1[1] = 0.0
Coordinates1[2] = 0.0
Radius1 = 0.51

Coordinates2 = KratosMultiphysics.Array3()
Coordinates2[0] = 1.0
Coordinates2[1] = 0.0
Coordinates2[2] = 0.0
Radius2 = 0.51


self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates1, Radius1, 1)
#Node1 = DEM1.GetNodes()[0]
self.ParticleCreatorDestructor.FEMDEM_CreateSphericParticle(Coordinates2, Radius2, 2)
#Node2 = DEM2.GetNodes()[0]

# Update the spheres model part
self.SpheresModelPart = self.ParticleCreatorDestructor.GetSpheresModelPart()
#counter = 0

for node in self.SpheresModelPart.Nodes:
	counter += 1
	#print(node.X)
	if (counter == 1):
		node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0.5)
	else:
		node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, -0.5)

#print(math.sqrt(25))
#print(self.SpheresModelPart.Nodes[0].X)
#print(self.SpheresModelPart.Nodes)

#Wait()
for element in self.SpheresModelPart.Elements:
	counter += 1
	print(element.Properties[KratosMultiphysics.YOUNG_MODULUS])

print("cont", counter)
Wait()

#Node1.SetSolutionStepValue(VELOCITY_X, 0.5)
#Node2.SetSolutionStepValue(VELOCITY_X, -0.5)


cont = 0.0
for node in self.SpheresModelPart.Nodes:
	cont += 1

print(cont)
Wait()


				for FEM_elem in self.FEM_Solution.main_model_part.Elements:
			print(FEM_elem.Id)
			for i in range(0,3):
				print(FEM_elem.GetNodes()[i].Id)
			Wait()
		
				is_active = True
		for FEM_elem in self.FEM_Solution.main_model_part.Elements:
			if FEM_elem.IsDefined(KratosMultiphysics.ACTIVE):
				is_active = FEM_elem.Is(KratosMultiphysics.ACTIVE)
				print(is_active)
				print(FEM_elem.Id)
		
		
		for FEM_elem in self.FEM_Solution.main_model_part.Elements:
			damage = FEM_elem.GetValue(KratosFemDem.DAMAGE_ELEMENT)
			if damage > 0.0:
				print(FEM_elem.Id)
				print(damage)
				print("*******")
		
		
		for node in self.FEM_Solution.main_model_part.Nodes:
			if node.Id == 3:
				#node.SetValue(KratosFemDem.IS_DEM, True)
				#node.SetValue(KratosFemDem.DEM_RADIUS, 1.22)
				print(self.GetNodeCoordinates(node))
				Wait()
		
		
		for node in self.FEM_Solution.main_model_part.Nodes:
			#print(node.GetValue(KratosFemDem.IS_DEM))
			#print(node.GetValue(KratosFemDem.DEM_RADIUS))
			#print(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X))
			print(node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
		Wait()
		

		
		for i in range(0,5):
			print("arriba", i)
			for j in range (0,5):
				print("debajo", j)
				if i == j:
					break
		

		# just for testing ->Remove
		#self.FEM_Solution.GraphicalOutputPrintOutput()
		# *********************** 

		#print(self.FEM_Solution.main_model_part.GetSubModelPart("Solid_Displacement-auto-1").Nodes
		#Wait()

		for node in self.FEM_Solution.main_model_part.GetSubModelPart("Solid_Displacement-auto-1").Nodes:
			print(node.Id)
		Wait()
		
		for node in self.FEM_Solution.main_model_part.Nodes:
			print(str(node.Id) + "  " + str(node.IsFixed(KratosMultiphysics.DISPLACEMENT_X)))
		Wait()
		# just for debug
		
		for node in self.FEM_Solution.main_model_part.Nodes:
			print(str(node.Id) + "  " + str(node.X) + "  " + str(node.Y))
		

		
		for elem in self.FEM_Solution.main_model_part.Elements:
			if elem.GetNodes()[0].Id == 55 or elem.GetNodes()[1].Id == 55 or elem.GetNodes()[2].Id == 55:
				print(str(elem.Id) + " " + str(elem.GetNodes()[0].Id)+ " " + str(elem.GetNodes()[1].Id)+ " " + str(elem.GetNodes()[2].Id))
				#Wait()
			else:
				print(str(elem.Id) + "  no node 55  " + str(elem.GetNodes()[0].Id)+ " " + str(elem.GetNodes()[1].Id)+ " " + str(elem.GetNodes()[2].Id))
		Wait() 
		# just for debug


				#Wait()
		'''print("nodos fijados...")
		for node in self.FEM_Solution.main_model_part.Nodes:
			if node.IsFixed(KratosMultiphysics.DISPLACEMENT_X) == True:
				print(node.Id, node.IsFixed(KratosMultiphysics.DISPLACEMENT_X))'''

		'''for node in self.FEM_Solution.main_model_part.GetSubModelPart("Solid_Displacement-auto-1").Nodes:
			print(node.Id)'''


		#print(self.FEM_Solution.main_model_part.GetSubModelPart("Solid_Displacement-auto-1"))
		'''for node in self.FEM_Solution.main_model_part.GetSubModelPart("Solid_Displacement-auto-1").Nodes:
			print(node.Id)
		Wait()'''