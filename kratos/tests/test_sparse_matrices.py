import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import numpy as np

try:
    import KratosMultiphysics.scipy_conversion_tools
    missing_scipy = False
except ImportError as e:
    missing_scipy = True

class TestSparseMatrixInterface(KratosUnittest.TestCase):

    def test_scalar_assembly(self):
        G = KratosMultiphysics.SparseContiguousRowGraph(3)
        G.AddEntry(0,0)
        G.AddEntry(2,2)
        G.AddEntry(1,2)
        G.Finalize()

        A = KratosMultiphysics.CsrMatrix(G)
        #5 0 0
        #0 0 7
        #0 0 6
        A.BeginAssemble()
        A.AssembleEntry(5.0,0,0)
        A.AssembleEntry(6.0,2,2)
        A.AssembleEntry(7.0,1,2)
        A.FinalizeAssemble()

        x = KratosMultiphysics.SystemVector(3)
        x.SetValue(1.0)
        y = KratosMultiphysics.SystemVector(3)
        y.SetValue(0.0)
        A.SpMV(x,y)
        self.assertEqual(y[0],5.0)
        self.assertEqual(y[1],7.0)
        self.assertEqual(y[2],6.0)

        y = KratosMultiphysics.SystemVector(3)
        y.SetValue(0.0)
        A.TransposeSpMV(x,y)
        self.assertEqual(y[0],5.0)
        self.assertEqual(y[1],0.0)
        self.assertEqual(y[2],13.0)

    def test_matrix_assembly(self):
        #data to be assembled
        values = np.array([[1.0,-1.0],[-1.0,1.0]])

        x = KratosMultiphysics.SystemVector(3)
        x[0] = 1.0
        x[1] = 2.0
        x[2] = 3.0
        y = KratosMultiphysics.SystemVector(3)
        y.SetValue(0.0)

        #construction of matrix graph (by vector input)
        G = KratosMultiphysics.SparseContiguousRowGraph(3)
        G.AddEntries(np.array([0,1])) #for square matrices
        G.AddEntries([1,2],[1,2]) #thos would allow non square input
        G.Finalize()

        #assembling matrix
        A = KratosMultiphysics.CsrMatrix(G)
        A.BeginAssemble()
        A.Assemble(values,[0,1]) #input by list
        A.Assemble(values,np.array([1,2])) #input by numpy array
        A.FinalizeAssemble()

        A.SpMV(x,y)

        validation_data = [ 2., -3.,  1., -3.,  6., -3.,  1., -3.,  2.]
        validation_index2 = [0, 1, 2, 0, 1, 2, 0, 1, 2]
        validation_index1 = [0, 3, 6, 9]

        B = A.SpMM(A)

        for i in range(len(validation_data)):
            self.assertEqual(B.value_data()[i], validation_data[i])
            self.assertEqual(B.index2_data()[i], validation_index2[i])

        for i in range(len(validation_index1)):
            self.assertEqual(B.index1_data()[i], validation_index1[i])

        B = A@A
        for i in range(len(validation_data)):
            self.assertEqual(B.value_data()[i], validation_data[i])
            self.assertEqual(B.index2_data()[i], validation_index2[i])

        for i in range(len(validation_index1)):
            self.assertEqual(B.index1_data()[i], validation_index1[i])

        # the following should be added back in case scipy support is enabled in testing
        if not missing_scipy:
             #test conversion to scipy matrix
            Ascipy = KratosMultiphysics.scipy_conversion_tools.to_csr(A)

            #verify compatibility with numpy arrays
            xscipy = np.array(x).T
            yscipy = Ascipy@x

            self.assertEqual(yscipy[0],-1.0)
            self.assertEqual(yscipy[1],0.0)
            self.assertEqual(yscipy[2],1.0)

            B_scipy = Ascipy@Ascipy
            B_scipy.sort_indices() #it is crucial that the index are sorted to do the following comparison

            for i in range(len(validation_data)):
                self.assertEqual(B_scipy.data[i], validation_data[i])
                self.assertEqual(B_scipy.indices[i], validation_index2[i])

            for i in range(len(validation_index1)):
                self.assertEqual(B_scipy.indptr[i], validation_index1[i])

        #test transpose and TransposeSpMV
        x.SetValue(1.0)
        y.SetValue(0.0)
        A.TransposeSpMV(x,y)

        At = A.Transpose()
        y2 =  KratosMultiphysics.SystemVector(y.size())
        y2.SetValue(0.0)
        At.SpMV(x,y2)
        for i in range(y.Size()):
            self.assertEqual(y[i], y2[i], 1e-14)

        #test @ operators
        y = A@x

        #test operations
        x.SetValue(1.0)
        y.SetValue(2.0)
        c = 2.0*x+y*2.0 - x
        for i in range(y.Size()):
            self.assertEqual(c[i], 5.0)

    def test_get_equation_id_csr_indices_1(self):
        # Create model part
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart("Main")
        mp.SetBufferSize(2)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        # Create nodes
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        mp.CreateNewNode(4, 1.0, 1.0, 0.0)

        # Create elements (Element2D3N)
        prop = mp.CreateNewProperties(1)
        mp.CreateNewElement("DistanceCalculationElementSimplex2D3N", 1, [1, 2, 3], prop)
        mp.CreateNewElement("DistanceCalculationElementSimplex2D3N", 2, [2, 4, 3], prop) # Mesh of 2 triangles

        # Add DOFs
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISTANCE)

        # Initialize DOFs (EquationIds)
        for node in mp.Nodes:
            node.GetDof(KratosMultiphysics.DISTANCE).EquationId = node.Id - 1

        # Build sparse matrix graph
        graph = KratosMultiphysics.SparseContiguousRowGraph(mp.NumberOfNodes())
        for elem in mp.Elements:
            equation_id_vector = elem.EquationIdVector(mp.ProcessInfo)
            graph.AddEntries(equation_id_vector, equation_id_vector)
        graph.Finalize()

        # Create Matrix
        A = KratosMultiphysics.CsrMatrix(graph)

        # Get the elemental contributions CSR indices
        elem_csr_indices = A.GetEquationIdCsrIndices(mp.Elements, mp.ProcessInfo)

        # Check results
        shape = elem_csr_indices.to_numpy().shape
        self.assertEqual(shape[0], 2)
        self.assertEqual(shape[1], 3)
        self.assertEqual(shape[2], 3)

        data = elem_csr_indices.to_numpy().flatten()
        expected_data = [0, 1, 2, 3, 4, 5, 7, 8, 9,
                         4, 6, 5, 11, 13, 12, 8, 10, 9]
        for i, j in zip(data, expected_data):
            self.assertEqual(i, j)

    def test_get_equation_id_csr_indices_2(self):
        # Create elements connectivities
        connectivities = KratosMultiphysics.IntNDData(np.array([[0,1,2],[1,3,2]]))

        # Build sparse matrix graph
        connectivities_data = connectivities.to_numpy()
        graph = KratosMultiphysics.SparseContiguousRowGraph(4)
        for i in range(connectivities_data.shape[0]):
            graph.AddEntries(connectivities_data[i], connectivities_data[i])
        graph.Finalize()

        # Create Matrix
        A = KratosMultiphysics.CsrMatrix(graph)

        # Get the elemental contributions CSR indices
        elem_csr_indices = A.GetEquationIdCsrIndices(connectivities)

        # Check results
        shape = elem_csr_indices.to_numpy().shape
        self.assertEqual(shape[0], 2)
        self.assertEqual(shape[1], 3)
        self.assertEqual(shape[2], 3)

        data = elem_csr_indices.to_numpy().flatten()
        expected_data = [0, 1, 2, 3, 4, 5, 7, 8, 9,
                         4, 6, 5, 11, 13, 12, 8, 10, 9]
        for i, j in zip(data, expected_data):
            self.assertEqual(i, j)

    def test_assemble_with_csr_indices(self):
        # Create model part
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart("Main")
        mp.SetBufferSize(2)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        # Create nodes
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        mp.CreateNewNode(4, 1.0, 1.0, 0.0)

        # Create elements (Element2D3N)
        prop = mp.CreateNewProperties(1)
        mp.CreateNewElement("DistanceCalculationElementSimplex2D3N", 1, [1, 2, 3], prop)
        mp.CreateNewElement("DistanceCalculationElementSimplex2D3N", 2, [2, 4, 3], prop) # Mesh of 2 triangles

        # Add DOFs
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISTANCE)

        # Initialize DOFs (EquationIds)
        for node in mp.Nodes:
            node.GetDof(KratosMultiphysics.DISTANCE).EquationId = node.Id - 1

        # Build sparse matrix graph
        graph = KratosMultiphysics.SparseContiguousRowGraph(mp.NumberOfNodes())
        for elem in mp.Elements:
            equation_id_vector = elem.EquationIdVector(mp.ProcessInfo)
            graph.AddEntries(equation_id_vector, equation_id_vector)
        graph.Finalize()

        # Create Matrix
        A = KratosMultiphysics.CsrMatrix(graph)

        # Get the elemental contributions CSR indices
        elem_csr_indices = A.GetEquationIdCsrIndices(mp.Elements, mp.ProcessInfo)

        # Calculate and store the LHS contributions
        elem_lhs_contributions = KratosMultiphysics.DoubleNDData([mp.NumberOfElements(), 3, 3])
        aux_data = elem_lhs_contributions.to_numpy()
        for elem in mp.Elements:
            lhs = KratosMultiphysics.Matrix()
            rhs = KratosMultiphysics.Vector()
            elem.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            for i in range(3):
                for j in range(3):
                    aux_data[(elem.Id-1)][i][j] = lhs[i, j]

        # Assemble the contributions
        A.AssembleWithCsrIndices(elem_lhs_contributions, elem_csr_indices)

        # Check results
        data = A.value_data()
        expected_data = [0.1, -0.05, -0.05, -0.05, 0.1, 0.0, -0.05, -0.05, 0.0, 0.1, -0.05, -0.05, -0.05, 0.1]
        self.assertVectorAlmostEqual(data, expected_data)

    def test_assemble_with_equation_ids(self):
        # Create model part
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart("Main")
        mp.SetBufferSize(2)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

        # Create nodes
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, 1.0, 0.0, 0.0)
        mp.CreateNewNode(3, 0.0, 1.0, 0.0)
        mp.CreateNewNode(4, 1.0, 1.0, 0.0)

        # Create elements (Element2D3N)
        prop = mp.CreateNewProperties(1)
        mp.CreateNewElement("DistanceCalculationElementSimplex2D3N", 1, [1, 2, 3], prop)
        mp.CreateNewElement("DistanceCalculationElementSimplex2D3N", 2, [2, 4, 3], prop) # Mesh of 2 triangles

        # Add DOFs
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.DISTANCE)

        # Initialize DOFs (EquationIds)
        for node in mp.Nodes:
            node.GetDof(KratosMultiphysics.DISTANCE).EquationId = node.Id - 1

        # Build sparse matrix graph
        graph = KratosMultiphysics.SparseContiguousRowGraph(mp.NumberOfNodes())
        for elem in mp.Elements:
            equation_id_vector = elem.EquationIdVector(mp.ProcessInfo)
            graph.AddEntries(equation_id_vector, equation_id_vector)
        graph.Finalize()

        # Create Matrix
        A = KratosMultiphysics.CsrMatrix(graph)

        # Get the elemental equation ids matrix
        equation_ids_data = np.empty((mp.NumberOfElements(), 3), dtype=np.int32)
        for elem in mp.Elements:
            equation_id_vector = elem.EquationIdVector(mp.ProcessInfo)
            equation_ids_data[elem.Id-1, :] = equation_id_vector
        elem_equation_ids = KratosMultiphysics.IntNDData(equation_ids_data)

        # Calculate and store the LHS contributions
        elem_lhs_contributions = KratosMultiphysics.DoubleNDData([mp.NumberOfElements(), 3, 3])
        aux_data = elem_lhs_contributions.to_numpy()
        for elem in mp.Elements:
            lhs = KratosMultiphysics.Matrix()
            rhs = KratosMultiphysics.Vector()
            elem.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
            for i in range(3):
                for j in range(3):
                    aux_data[(elem.Id-1)][i][j] = lhs[i, j]

        # Assemble the contributions
        A.AssembleWithEquationIds(elem_lhs_contributions, elem_equation_ids)

        # Check results
        data = A.value_data()
        expected_data = [0.1, -0.05, -0.05, -0.05, 0.1, 0.0, -0.05, -0.05, 0.0, 0.1, -0.05, -0.05, -0.05, 0.1]
        self.assertVectorAlmostEqual(data, expected_data)

if __name__ == '__main__':
    KratosUnittest.main()
