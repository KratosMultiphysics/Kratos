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
    
class ShellQ4ThickBendingRollUpTests(KratosUnittest.TestCase):

    def setUp(self):
      self.file_name = "Shell_Q4_Thick__BendingRollUp_test/Shell_Q4_Thick__BendingRollUp_test"
      # Initialize GiD  I/O
      input_file_name = GetFilePath(self.file_name)
      parameter_file = open(self.file_name +"_parameters.json",'r')
      ProjectParameters = Parameters( parameter_file.read())
      # Creating the model part
      self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_ShellQ4ThickBendingRollUpTests(self):
        self.test.Solve()

    def tearDown(self):
        pass
    
class ShellQ4ThickDrillingRollUpTests(KratosUnittest.TestCase):

    def setUp(self):
      self.file_name = "Shell_Q4_Thick__DrillingRollUp_test/Shell_Q4_Thick__DrillingRollUp_test"
      # Initialize GiD  I/O
      input_file_name = GetFilePath(self.file_name)
      parameter_file = open(self.file_name +"_parameters.json",'r')
      ProjectParameters = Parameters( parameter_file.read())
      # Creating the model part
      self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_ShellQ4ThickDrillingRollUpTests(self):
        self.test.Solve()

    def tearDown(self):
        pass
    
    
class ShellT3IsotropicScordelisTests(KratosUnittest.TestCase):

    def setUp(self):
      self.file_name = "Shell_T3_Isotropic_Scordelis_test/Shell_T3_Isotropic_Scordelis_test"
      # Initialize GiD  I/O
      input_file_name = GetFilePath(self.file_name)
      parameter_file = open(self.file_name +"_parameters.json",'r')
      ProjectParameters = Parameters( parameter_file.read())
      # Creating the model part
      self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_ShellT3IsotropicScordelisTests(self):
        self.test.Solve()

    def tearDown(self):
        pass
    
    
class ShellT3ThinBendingRollUpTests(KratosUnittest.TestCase):

    def setUp(self):
      self.file_name = "Shell_T3_Thin__BendingRollUp_test/Shell_T3_Thin__BendingRollUp_test"
      # Initialize GiD  I/O
      input_file_name = GetFilePath(self.file_name)
      parameter_file = open(self.file_name +"_parameters.json",'r')
      ProjectParameters = Parameters( parameter_file.read())
      # Creating the model part
      self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_ShellT3ThinBendingRollUpTests(self):
        self.test.Solve()

    def tearDown(self):
        pass


class ShellT3ThinDrillingRollUpTests(KratosUnittest.TestCase):

    def setUp(self):
      self.file_name = "Shell_T3_Thin__DrillingRollUp_test/Shell_T3_Thin__DrillingRollUp_test"
      # Initialize GiD  I/O
      input_file_name = GetFilePath(self.file_name)
      parameter_file = open(self.file_name +"_parameters.json",'r')
      ProjectParameters = Parameters( parameter_file.read())
      # Creating the model part
      self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_ShellT3ThinDrillingRollUpTests(self):
        self.test.Solve()

    def tearDown(self):
        pass
