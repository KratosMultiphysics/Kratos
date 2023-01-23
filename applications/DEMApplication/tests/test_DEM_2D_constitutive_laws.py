import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import KratosMultiphysics.kratos_utilities as kratos_utils
import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class DEM2DConstitutiveLawsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_constitutive_laws_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def InitializeSolutionStep(self):

        if self.time > self.time_to_heal_and_set_indentation and not self.already_called_heal_and_setting_indentation:
            self._GetSolver().cplusplus_strategy.HealAllBonds()
            KratosMultiphysics.DEMApplication.ParallelBondUtilities().SetCurrentIndentationAsAReferenceInParallelBonds(self.spheres_model_part)
            self.already_called_heal_and_setting_indentation = True

        super().InitializeSolutionStep()

    def Initialize(self):
        super().Initialize()
        self.print_pdf_files = False
        self.already_called_heal_and_setting_indentation = False
        self.time_to_heal_and_set_indentation = 2e-7
        self._GetSolver().cplusplus_strategy.BreakAllBonds()

    def Finalize(self):
        if self.print_pdf_files:
            self.PrintDebugGraphs()
        super().Finalize()

    def PrintDebugGraphs(self):

        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf
        pdf_name = self.GetProblemTypeFileName() + '.pdf'

        pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name)
        small = 5
        plt.rc('figure', titlesize=small); plt.rc('font', size=small); plt.rc('axes', titlesize=small); plt.rc('axes', labelsize=small)
        plt.rc('xtick', labelsize=small); plt.rc('ytick', labelsize=small); plt.rc('legend', fontsize=small)
        wspace = 0.3
        hspace = 0.3

        #Normal forces
        X, Y1, Y2, Y3, Y4, Y5, Y6 = [], [], [], [], [], [], []
        Y7, Y8, Y9, Y10, Y11, Y12 = [], [], [], [], [], []
        Y13, Y14, Y15, Y16, Y17, Y18 = [], [], [], [], [], []
        Y19, Y20, Y21, Y22 = [], [], [], []
        with open(os.path.join(self.GetMainPath(), 'nl.txt'), 'r') as normal:
            for line in normal:
                values = [float(s) for s in line.split()]
                X.append(values[0])
                Y1.append(values[1]); Y2.append(values[2]); Y3.append(values[3]); Y4.append(values[4]); Y5.append(values[5]); Y6.append(values[6])
                Y7.append(values[7]); Y8.append(values[8]); Y9.append(values[9]); Y10.append(values[10]); Y11.append(values[11]); Y12.append(values[12])
                Y13.append(values[13]); Y14.append(values[14]); Y15.append(values[15]); Y16.append(values[16]); Y17.append(values[17]); Y18.append(values[18])
                Y19.append(values[19]); Y20.append(values[20]); Y21.append(values[21]); Y22.append(values[22])
        plt.figure(1)
        plt.subplot(3, 2, 1); plt.plot(X, Y1, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) indentation'); plt.grid()
        plt.subplot(3, 2, 2); plt.plot(X, Y2, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) LocalElasticContactForce[2]'); plt.grid()
        plt.subplot(3, 2, 3); plt.plot(X, Y3, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) limit_force'); plt.grid()
        plt.subplot(3, 2, 4); plt.plot(X, Y4, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) delta_accumulated'); plt.grid()
        plt.subplot(3, 2, 5); plt.plot(X, Y5, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) returned_by_mapping_force'); plt.grid()
        plt.subplot(3, 2, 6); plt.plot(X, Y6, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) kn_updated'); plt.grid()
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        pdf.savefig(1)
        plt.figure(2)
        plt.subplot(3, 2, 1); plt.plot(X, Y7, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) mDamageNormal'); plt.grid()
        plt.subplot(3, 2, 2); plt.plot(X, Y8, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) failure_type'); plt.grid()
        plt.subplot(3, 2, 3); plt.plot(X, Y9, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) current_normal_force_module'); plt.grid()
        plt.subplot(3, 2, 4); plt.plot(X, Y10, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) mDamageTangential'); plt.grid()
        plt.subplot(3, 2, 5); plt.plot(X, Y11, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) BondedLocalElasticContactForce2'); plt.grid()
        plt.subplot(3, 2, 6); plt.plot(X, Y12, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) mUnbondedLocalElasticContactForce2'); plt.grid()
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        pdf.savefig(2)
        plt.figure(3)
        plt.subplot(3, 2, 1); plt.plot(X, Y13, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) kn_el'); plt.grid()
        plt.subplot(3, 2, 2); plt.plot(X, Y14, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) mDamageEnergyCoeff'); plt.grid()
        plt.subplot(3, 2, 3); plt.plot(X, Y15, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) initial_limit_force'); plt.grid()
        plt.subplot(3, 2, 4); plt.plot(X, Y16, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) mUnbondedNormalElasticConstant'); plt.grid()
        plt.subplot(3, 2, 5); plt.plot(X, Y17, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) LocalElasticContactTension'); plt.grid()
        plt.subplot(3, 2, 6); plt.plot(X, Y18, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) limit_tension'); plt.grid()
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        pdf.savefig(3)
        plt.figure(4)
        plt.subplot(2, 2, 1); plt.plot(X, Y19, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) returned_by_mapping_tension'); plt.grid()
        plt.subplot(2, 2, 2); plt.plot(X, Y20, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) current_normal_tension_module'); plt.grid()
        plt.subplot(2, 2, 3); plt.plot(X, Y21, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) BondedLocalElasticContactTension2'); plt.grid()
        plt.subplot(2, 2, 4); plt.plot(X, Y22, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(N) bonded_indentation'); plt.grid()
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        pdf.savefig(4)

        #Tangential forces
        X, Y1, Y2, Y3, Y4, Y5, Y6 = [], [], [], [], [], [], []
        Y7, Y8, Y9, Y10, Y11, Y12     = [], [], [], [], [], []
        Y13, Y14, Y15, Y16, Y17, Y18  = [], [], [], [], [], []
        Y19, Y20, Y21, Y22, Y23, Y24  = [], [], [], [], [], []
        Y25, Y26, Y27, Y33, Y34, Y35  = [], [], [], [], [], []
        with open(os.path.join(self.GetMainPath(), 'tg.txt'), 'r') as tangential:
            for line in tangential:
                values = [float(s) for s in line.split()]
                X.append(values[0])
                Y1.append(values[1]); Y2.append(values[2]); Y3.append(values[3]); Y4.append(values[4]); Y5.append(values[5]); Y6.append(values[6])
                Y7.append(values[7]); Y8.append(values[8]); Y9.append(values[9]); Y10.append(values[10]); Y11.append(values[11]); Y12.append(values[12])
                Y13.append(values[13]); Y14.append(values[14]); Y15.append(values[15]); Y16.append(values[16]); Y17.append(values[17]); Y18.append(values[18])
                Y19.append(values[19]); Y20.append(values[20]); Y21.append(values[21]); Y22.append(values[22]); Y23.append(values[23]); Y24.append(values[24])
                Y25.append(values[25]); Y26.append(values[26]); Y27.append(values[27]); Y33.append(values[33]); Y34.append(values[34]); Y35.append(values[35])
        plt.figure(5)
        plt.subplot(3, 2, 1); plt.plot(X, Y1, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) failure_type'); plt.grid()
        plt.subplot(3, 2, 2); plt.plot(X, Y2, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) tau_strength'); plt.grid()
        plt.subplot(3, 2, 3); plt.plot(X, Y3, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) kt_updated'); plt.grid()
        plt.subplot(3, 2, 4); plt.plot(X, Y4, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) sliding'); plt.grid()
        plt.subplot(3, 2, 5); plt.plot(X, Y5, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) contact_sigma'); plt.grid()
        plt.subplot(3, 2, 6); plt.plot(X, Y6, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) mDamageNormal'); plt.grid()
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        pdf.savefig(5)
        plt.figure(6)
        plt.subplot(3, 2, 1); plt.plot(X, Y7, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) contact_tau'); plt.grid()
        plt.subplot(3, 2, 2); plt.plot(X, Y8, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) max_admissible_shear_force'); plt.grid()
        plt.subplot(3, 2, 3); plt.plot(X, Y9, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) mDamageTangential'); plt.grid()
        plt.subplot(3, 2, 4); plt.plot(X, Y10, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) LocalElasticContactForce[0]'); plt.grid()
        plt.subplot(3, 2, 5); plt.plot(X, Y11, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) OldBondedLocalElasticContactForce[0]'); plt.grid()
        plt.subplot(3, 2, 6); plt.plot(X, Y12, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) mUnbondedLocalElasticContactForce2'); plt.grid()
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        pdf.savefig(6)
        plt.figure(7)
        plt.subplot(3, 2, 1); plt.plot(X, Y13, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) BondedLocalElasticContactForce[0]'); plt.grid()
        plt.subplot(3, 2, 2); plt.plot(X, Y14, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) mBondedScalingFactor'); plt.grid()
        plt.subplot(3, 2, 3); plt.plot(X, Y15, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) UnbondedLocalElasticContactForce[0]'); plt.grid()
        plt.subplot(3, 2, 4); plt.plot(X, Y16, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) LocalDeltDisp[0]'); plt.grid()
        plt.subplot(3, 2, 5); plt.plot(X, Y17, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) mUnbondedScalingFactor'); plt.grid()
        plt.subplot(3, 2, 6); plt.plot(X, Y18, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) OldLocalElasticContactForce[0]'); plt.grid()
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        pdf.savefig(7)
        plt.figure(8)
        plt.subplot(3, 2, 1); plt.plot(X, Y19, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) LocalDeltDisp[1]'); plt.grid()
        plt.subplot(3, 2, 2); plt.plot(X, Y20, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) LocalElasticContactForce[1]'); plt.grid()
        plt.subplot(3, 2, 3); plt.plot(X, Y21, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) BondedLocalElasticContactForce[1]'); plt.grid()
        plt.subplot(3, 2, 4); plt.plot(X, Y22, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) UnbondedLocalElasticContactForce[1]'); plt.grid()
        plt.subplot(3, 2, 5); plt.plot(X, Y23, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) OldBondedLocalElasticContactForce[1]'); plt.grid()
        plt.subplot(3, 2, 6); plt.plot(X, Y24, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) OldUnbondedLocalElasticContactForce[1]'); plt.grid()
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        pdf.savefig(8)
        plt.figure(9)
        plt.subplot(3, 2, 1); plt.plot(X, Y25, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) delta_accumulated'); plt.grid()
        plt.subplot(3, 2, 2); plt.plot(X, Y26, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) current_tangential_force_module'); plt.grid()
        plt.subplot(3, 2, 3); plt.plot(X, Y27, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) returned_by_mapping_force'); plt.grid()
        plt.subplot(3, 2, 4); plt.plot(X, Y33, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) mDamageEnergyCoeff'); plt.grid()
        plt.subplot(3, 2, 5); plt.plot(X, Y34, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) LocalElasticContactForce[2]'); plt.grid()
        plt.subplot(3, 2, 6); plt.plot(X, Y35, '-', color='blue'); plt.xlabel('time (s)'); plt.ylabel('(T) indentation'); plt.grid()
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        pdf.savefig(9)
        pdf.close()
        plt.close('all')

class DEM2DConstitutiveLaws(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_DEM2D_ConstitutiveLaws1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_constitutive_laws_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM1.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2DConstitutiveLawsTestSolution, model, parameters_file_name, 1)

    @classmethod
    def test_DEM2D_ConstitutiveLaws2(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_constitutive_laws_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM2.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2DConstitutiveLawsTestSolution, model, parameters_file_name, 1)

    @classmethod
    def test_DEM2D_ConstitutiveLaws3(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "DEM2D_constitutive_laws_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM3.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(DEM2DConstitutiveLawsTestSolution, model, parameters_file_name, 1)

    def tearDown(self):
        file_to_remove = os.path.join("DEM2D_constitutive_laws_tests_files", "TimesPartialRelease")
        kratos_utils.DeleteFileIfExisting(GetFilePath(file_to_remove))
        file_to_remove = os.path.join(this_working_dir_backup, "nl.txt")
        os.remove(file_to_remove)
        file_to_remove = os.path.join(this_working_dir_backup, "tg.txt")
        os.remove(file_to_remove)
        os.chdir(this_working_dir_backup)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
