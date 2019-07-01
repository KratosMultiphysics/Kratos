import time
import KratosMultiphysics as km
import KratosMultiphysics.FluidDynamicsApplication as kfd
import KratosMultiphysics.ChimeraApplication as kchim
try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError:
    have_external_solvers = False

from fluid_chimera_analysis import FluidChimeraAnalysis
from chimera_analysis_base_test import ChimeraAnalysisBaseTest

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities

@UnitTest.skipUnless(have_external_solvers,"Missing required application: ExternalSolversApplication")
class FlowOverCylinderMonolithic(ChimeraAnalysisBaseTest):
    def test_MonolithicFlowOverCylinder(self):
        start = time.clock()
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test

        work_folder = "flow_over_cylinder_monolithic"
        settings_file_name = "flow_over_cylinder_monolithic.json"
        with UnitTest.WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
        end = time.clock()
        print("Time taken for monolithic chimera simulation",end-start)

class FlowOverCylinderFractionalStep(ChimeraAnalysisBaseTest):
    def test_FractionalStepFlowOverCylinder(self):
        start = time.clock()
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        work_folder = "flow_over_cylinder_fractionalstep"
        settings_file_name = "flow_over_cylinder_fractionalstep.json"
        with UnitTest.WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
        end = time.clock()
        print("Time taken for fractional step chimera simulation",end-start)

class MonolithicMultiPatch(ChimeraAnalysisBaseTest):
    def test_MultipleOverlappingPatchMonolithic(self):
        start = time.clock()
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        work_folder = "chimera_multiple_overlapping_patch"
        settings_file_name = "test_chimera_multiple_overlapping_simple_ProjectParameters.json"
        with UnitTest.WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
        end = time.clock()
        print("Time taken for Multiple overlapping chimera simulation using Monolithic solver ",end-start)

class FractionalStepMultiPatch(ChimeraAnalysisBaseTest):
    def test_MultipleOverlappingPatchFractionalStep(self): #TODO: Check and correct this
        start = time.clock()
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        work_folder = "chimera_multiple_overlapping_patch"
        settings_file_name = "test_chimera_multiple_overlapping_simple_ProjectParameters.json"
        with UnitTest.WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
        end = time.clock()
        print("Time taken for Multiple overlapping chimera simulation using fractional step solver ",end-start)


class FlowOverCrossMonolithic(ChimeraAnalysisBaseTest):
    def test_FlowOverCrossMonolithic(self):
        start = time.clock()
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        work_folder = "flow_over_cross_monolithic"
        settings_file_name = "test_chimera_multiple_overlapping_simple_ProjectParameters.json"
        with UnitTest.WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
        end = time.clock()
        print("Time taken for flow over cross using monolithic solver ",end-start)

class FlowOverCrossFractionalStep(ChimeraAnalysisBaseTest):
    def test_FlowOverCrossFractionalStep(self):
        start = time.clock()
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        work_folder = "flow_over_cross_fractional_step"
        settings_file_name = "test_chimera_multiple_overlapping_simple_ProjectParameters.json"
        with UnitTest.WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
        end = time.clock()
        print("Time taken for flow over cross using fractionalstep solver ",end-start)