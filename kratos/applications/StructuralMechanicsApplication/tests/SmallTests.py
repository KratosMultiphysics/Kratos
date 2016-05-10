# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

import Kratos_Execute_Solid_Test as Execute_Test

def GetFilePath(fileName):
  return os.path.dirname(__file__) + "/" + fileName
  
class DynamicBossakTests(KratosUnittest.TestCase):

    def setUp(self):
      self.file_name = "dynamic_test/dynamic_bossak_test"
      # Initialize GiD  I/O
      input_file_name = GetFilePath(self.file_name)
      parameter_file = open(self.file_name +"_parameters.json",'r')
      ProjectParameters = Parameters( parameter_file.read())
      # Creating the model part
      self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)
      
    def test_Bossak(self):
      self.test.Solve()
      
    def tearDown(self):
        pass
    
class DynamicNewmarkTests(KratosUnittest.TestCase):

    def setUp(self):
      self.file_name = "dynamic_test/dynamic_newmark_test"
      # Initialize GiD  I/O
      input_file_name = GetFilePath(self.file_name)
      parameter_file = open(self.file_name +"_parameters.json",'r')
      ProjectParameters = Parameters( parameter_file.read())
      # Creating the model part
      self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)
      
    def test_Newmark(self):
      self.test.Solve()
      
    def tearDown(self):
        pass
    
class SprismMembranePatchTests(KratosUnittest.TestCase):

    def setUp(self):
      self.file_name = "sprism_test/patch_membrane_test"
      # Initialize GiD  I/O
      input_file_name = GetFilePath(self.file_name)
      parameter_file = open(self.file_name +"_parameters.json",'r')
      ProjectParameters = Parameters( parameter_file.read())
      # Creating the model part
      self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_MembranePatch(self):
        self.test.Solve()

    def tearDown(self):
        pass
    
class SprismBendingPatchTests(KratosUnittest.TestCase):

    def setUp(self):
      self.file_name = "sprism_test/patch_bending_test"
      # Initialize GiD  I/O
      input_file_name = GetFilePath(self.file_name)
      parameter_file = open(self.file_name +"_parameters.json",'r')
      ProjectParameters = Parameters( parameter_file.read())
      # Creating the model part
      self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_BendingPatch(self):
        self.test.Solve()

    def tearDown(self):
        pass