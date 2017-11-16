import os, random, h5py
from KratosMultiphysics import *
from KratosMultiphysics.HDF5Application import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import xdmf_utils

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class TestCase(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def _remove_file(self, file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)
    
    def _initialize_model_part(self, model_part):
        num_tri_elems = 10
        num_quad_elems = 15
        num_tri_conds = 5
        num_quad_conds = 10
        num_nodes = (num_tri_elems + 2) + (num_quad_elems + 2)
        # Add variables.
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        # Create nodes.
        for i in range(num_nodes):
            x = 0.0
            y = 0.0
            z = 0.0
            model_part.CreateNewNode(i + 1, x, y, z)
        # Create triangle elements.
        for i in range(num_tri_elems):
            elem_id = i + 1
            node_ids = [i + 1, i + 2, i + 3]
            prop_id = 1
            prop = model_part.GetProperties()[prop_id]
            model_part.CreateNewElement("Element2D3N", elem_id, node_ids, prop)
        # Create quad elements.
        for i in range(num_quad_elems):
            elem_id = num_tri_elems + i + 1
            node_ids = [num_tri_elems + i + 1, num_tri_elems + i + 2, num_tri_elems + i + 3, num_tri_elems + i + 4]
            prop_id = 1
            prop = model_part.GetProperties()[prop_id]
            model_part.CreateNewElement("Element2D4N", elem_id, node_ids, prop)
        # Create triangle conditions.
        for i in range(num_tri_conds):
            cond_id = i + 1
            node_ids = [i + 1, i + 2, i + 3]
            prop_id = 1
            prop = model_part.GetProperties()[prop_id]
            model_part.CreateNewCondition("SurfaceCondition3D3N", cond_id, node_ids, prop)
        # Create quad conditions.
        for i in range(num_quad_conds):
            cond_id = num_tri_conds + i + 1
            node_ids = [num_tri_conds + i + 1, num_tri_conds + i + 2, num_tri_conds + i + 3, num_tri_conds + i + 4]
            prop_id = 1
            prop = model_part.GetProperties()[prop_id]
            model_part.CreateNewCondition("SurfaceCondition3D4N", cond_id, node_ids, prop)
        model_part.SetBufferSize(2)

    def _get_file(self):
        params = Parameters("""
        {
            "file_name" : "test_xdmf_output.h5",
            "file_access_mode" : "truncate",
            "file_driver" : "sec2"
        }""")
        return HDF5FileSerial(params)

    def _get_model_part_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ModelData",
            "list_of_elements" : ["Element2D3N", "Element2D4N"],
            "list_of_conditions" : ["SurfaceCondition3D3N", "SurfaceCondition3D4N"]
        }""")
        return HDF5ModelPartIO(params, hdf5_file)

    def _get_nodal_solution_step_data_io(self, hdf5_file):
        params = Parameters("""
        {
            "partitioned" : false,
            "prefix" : "/ResultsData",
            "list_of_variables" : ["DISPLACEMENT", "DENSITY", "PARTITION_INDEX"]
        }""")
        return HDF5NodalSolutionStepDataIO(params, hdf5_file)

    def test_XdmfTopology(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model_part = ModelPart("test")
            self._initialize_model_part(model_part)
            kratos_hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(kratos_hdf5_file)
            hdf5_model_part_io.WriteModelPart(model_part)
            del hdf5_model_part_io, kratos_hdf5_file
            h5py_file = h5py.File("test_xdmf_output.h5", "r")
            xdmf_topology = xdmf_utils.XdmfTopology(h5py_file, '/ModelData/Elements/Element2D3N')
            xdmf_topology.Build()
            #xml_string = xdmf_utils.ET.tostring(xdmf_topology.GetXmlElement())
            self._remove_file("test_xdmf_output.h5")

    def test_XdmfGeometry(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model_part = ModelPart("test")
            self._initialize_model_part(model_part)
            kratos_hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(kratos_hdf5_file)
            hdf5_model_part_io.WriteModelPart(model_part)
            del hdf5_model_part_io, kratos_hdf5_file
            xdmf_utils.WriteXdmfParametricCoordinates("test_xdmf_output.h5", "/ModelData/Nodes/Local")
            h5py_file = h5py.File("test_xdmf_output.h5", "r")
            xdmf_geometry = xdmf_utils.XdmfGeometry(h5py_file, '/ModelData/Nodes/Local')
            xdmf_geometry.Build()
            self._remove_file("test_xdmf_output.h5")

    def test_XdmfNodalAttribute(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model_part = ModelPart("test")
            self._initialize_model_part(model_part)
            kratos_hdf5_file = self._get_file()
            hdf5_nodal_solution_step_data_io = self._get_nodal_solution_step_data_io(kratos_hdf5_file)
            hdf5_nodal_solution_step_data_io.WriteNodalResults(model_part.Nodes, 0)
            del hdf5_nodal_solution_step_data_io, kratos_hdf5_file
            h5py_file = h5py.File("test_xdmf_output.h5", "r")
            xdmf_scalar_attribute = xdmf_utils.XdmfNodalAttribute(h5py_file, '/ResultsData/NodalResults/DENSITY')
            xdmf_scalar_attribute.Build()
            #print(xdmf_utils.ET.tostring(xdmf_scalar_attribute.GetXmlElement()))
            xdmf_vector_attribute = xdmf_utils.XdmfNodalAttribute(h5py_file, '/ResultsData/NodalResults/DISPLACEMENT')
            xdmf_vector_attribute.Build()
            #print(xdmf_utils.ET.tostring(xdmf_vector_attribute.GetXmlElement()))
            self._remove_file("test_xdmf_output.h5")
    
    def test_XdmfModelPart(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model_part = ModelPart("test")
            self._initialize_model_part(model_part)
            kratos_hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(kratos_hdf5_file)
            hdf5_model_part_io.WriteModelPart(model_part)
            del hdf5_model_part_io, kratos_hdf5_file
            xdmf_utils.WriteXdmfParametricCoordinates("test_xdmf_output.h5", "/ModelData/Nodes/Local")
            h5py_file = h5py.File("test_xdmf_output.h5", "r")
            xdmf_model_part = xdmf_utils.XdmfModelPart(h5py_file, '/ModelData')
            xdmf_model_part.Build()
            #xdmf_utils.ET.ElementTree(xdmf_model_part.GetXmlElement()).write("test.xml")
            self._remove_file("test_xdmf_output.h5")

    def test_XdmfStaticResults(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model_part = ModelPart("test")
            self._initialize_model_part(model_part)
            kratos_hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(kratos_hdf5_file)
            hdf5_model_part_io.WriteModelPart(model_part)
            hdf5_nodal_solution_step_data_io = self._get_nodal_solution_step_data_io(kratos_hdf5_file)
            hdf5_nodal_solution_step_data_io.WriteNodalResults(model_part.Nodes, 0)
            del hdf5_nodal_solution_step_data_io, hdf5_model_part_io, kratos_hdf5_file
            xdmf_utils.WriteXdmfParametricCoordinates("test_xdmf_output.h5", "/ModelData/Nodes/Local")
            xdmf_results = xdmf_utils.XdmfStaticResults("test_xdmf_output.h5", "/ModelData", "test_xdmf_output.h5", "ResultsData")
            xdmf_results.Build()
            #xdmf_utils.ET.ElementTree(xdmf_results.GetXmlElement()).write("test.xml")
            self._remove_file("test_xdmf_output.h5")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
