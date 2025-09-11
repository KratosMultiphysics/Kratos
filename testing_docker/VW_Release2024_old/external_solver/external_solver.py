
"""
Author: David Schmölz
"""

import os
import KratosMultiphysics as Kratos


class ExternalSolver(object):

    def __init__(self, parameter_file_name=None):
        """
        Konstruktor des Externen Solver Interfaces.
        In diesem muss der Name des Solvers ('self.solver_name')
        definiert werden. Optional kann eine Parameterdatei
        mit zusätzlichen Einstellung für diesen Solver
        eingelesen werden.

        Z.b. als json Datei:
        parameter_file = open(parameter_file_name)
        self.parameter = json.load(parameter_file)

        Args:
            self
            parameter_file_name:    Optionale Parameterdatei (z.B. string)

        Returns:
            None
        """
        self.solver_name = None

    def GetCaseName(self) -> str:
        """
        Gibt den Namen des Solvers aus. Dieser wird verwendet
        um einen einzigartigen Pfad zu kreiieren ("{solver_counter}_{case_name}_iterations"),
        welcher dann Ordner für jeden einzelnen Optimierungsschritt enthält
        ('current_iteration_path'), in welchem die Solverdaten gespeichert sind.

        Z.b.:
        1_dummy_case_iterations         => Order für diesen Solvern
        |___1                           => Ordner Optimierungsschritt 1
        |   |   dummy_solver_file.fem   => Solver Datei
        |
        |___2                           => Ordner Optimierungsschritt 2
            |   dummy_solver_file.fem   => Solver Datei

        Args:
            self

        Returns:
            string: name of the solver case
        """
        NotImplementedError("ExternalSolver:: GetCaseName has to be implemented!")

    def GetOriginalFilePath(self) -> str:
        """
        Gibt den Namen des Ordners mit den Originaldateien des Solvers aus.

        Args:
            self

        Returns:
            string: name of original solver file path
        """
        NotImplementedError("ExternalSolver:: GetOriginalFilePath has to be implemented!")

    def InitializeBeforeOptimization(self) -> None:
        """
        In dieser optionalen Methode können Änderungen oder Initialisierungen des Solver
        (oder eventuell der Solverdateien) bevor die Kratos Optimierung startet
        vorgenommen werden.

        Args:
            self

        Returns:
            None
        """
        NotImplementedError("ExternalSolver:: InitializeBeforeOptimization has to be implemented!")

    def TranslateToKratosModel(self) -> None:
        """
        Übersetzen des Solver FE Modells in ein Kratos Modell.
        Das übersetzte Kratos Modell muss als .mdpa Datei
        ausgeschrieben werden.

        Z.b. werden Optistruct Elemente zu Kratos Conditions wie folgt
        übersetzt:
        {
            "CTRIA3" : "SurfaceCondition3D3N",
            "CQUAD4" : "SurfaceCondition3D4N",
            "CTRIA6" : "SurfaceCondition3D6N",
            "CQUAD8" : "SurfaceCondition3D8N"
        }

        Args:
            self

        Returns:
            None
        """
        NotImplementedError("ExternalSolver:: TranslateToKratosModel has to be implemented!")

    def CopyInputFilesToCurrentIteration(self,
                                         current_iteration_path: os.PathLike,
                                         previous_iteration_path: os.PathLike,
                                         original_path: os.PathLike,
                                         iteration: int) -> None:
        """
        Kopieren der Solver Dateien zu derzeitigen Iterationspfad ('current_iteration_path').
        Dabei können die Solver Dateien aus folgenden Pfaden verwendet werden:
        - vorherige Optimierungsiteration ('previous_iteration_path')
        - originale Dateien ('original_path')

        Args:
            current_iteration_path:     Pfad des derzeitigen Iterationsschrittes
            previous_iteration_path:    Pfad des vorherigen Iterationsschrittes
            original_path:              Pfad zu den originale Solver Dateien
            iteration:                  Zahl des derzeitigen Iterationsschrittes (erste Iteration = 1)

        Returns:
            None
        """
        NotImplementedError("ExternalSolver:: CopyInputFilesToCurrentIteration has to be implemented!")

    def UpdateMesh(self,
                   current_iteration_path: os.PathLike,
                   previous_iteration_path: os.PathLike,
                   original_path: os.PathLike,
                   iteration: int,
                   shape_update: dict,
                   thickness_update: dict,
                   shape: dict,
                   thickness: dict,
                   current_design: Kratos.Model()) -> None:
        """
        Aktualisieren des Solver FE Modells im derzeitigen
        Optimierungsiterationspfad ('current_iteration_path'),
        in welchem die Dateien im vorherigen Schritt
        ('CopyInputFilesToCurrentIteration') kopiert wurden.

        Folgende Daten werden von Kratos übergeben:

        - shape_update:     python Dictionary
                            Schlüssel:  Knoten Id
                            Werte:      Formänderung in globalen kartesischen Koordinaten
                            {Knoten Id: [Update X, Update Y, Update Z]}
                            {int: [float, float, float]}

        - thickness_update: python Dictionary
                            Schlüssel:  Element Id
                            Werte:      Schalendickenänderung
                            {Element Id: Update T}
                            {int: float}

        - shape:            python Dictionary
                            Schlüssel:  Knoten Id
                            Werte:      Knotenposition in globalen kartesischen Koordinaten
                            {Knoten Id: [X, Y, Z]}
                            {int: [float, float, float]}

        - thickness:        python Dictionary
                            Schlüssel:  Element Id
                            Werte:      Schalendicken
                            {Element Id: T}
                            {int: float}

        - current_design:   Kratos Optimierungs-ModelPart
                            Enthält alle Kratos Daten
                            z.B. Element/Conditions, Nodes, NodalVariables, etc.

        Args:
            current_iteration_path:     Pfad des derzeitigen Iterationsschrittes
            previous_iteration_path:    Pfad des vorherigen Iterationsschrittes
            original_path:              Pfad zu den originale Solver Dateien
            iteration:                  Zahl des derzeitigen Iterationsschrittes (erste Iteration = 1)
            shape_update:               Formänderung aus Kratos
            thickness_update:           Schalendickenänderung aus Kratos
            shape:                      Form aus Kratos
            thickness:                  Schalendicken aus Kratos
            current_design:             Derzeitiges Kratos Model

        Returns:
            None
        """
        NotImplementedError("ExternalSolver:: UpdateMesh has to be implemented!")

    def RunAnalysis(self,
                    current_iteration_path: os.PathLike,
                    current_design: Kratos.Model(),
                    iteration: int) -> None:
        """
        Methode in der die primäre und eventuell adjungierte Analyse
        des externen Solvers im derzeitigen Optimierungspfad laufen muss.
        Falls dieser Prozess fehlschlägt muss Kratos mit der Methode
        'self._StopKratosOptimizationRun()' gestoppt werden.
        Die Kratos Formoptimierung erwartet für jede Antwortfunktion, welche
        in der .json Datei der Kratos Optimierungsparameter spezifiert wurde
        den Wert der Funktion und Sensitivitäten dazu.
        Der Kratos Optimierungs-ModelPart 'current_design' kann optional
        verwendet werden.

        Args:
            current_iteration_path:     Pfad des derzeitigen Iterationsschrittes
            current_design:             Derzeitiges Kratos Model
            iteration:                  Zahl des derzeitigen Iterationsschrittes (erste Iteration = 1)

        Returns:
            None
        """
        NotImplementedError("ExternalSolver:: RunAnalysis has to be implemented!")

    def _StopKratosOptimizationRun() -> None:
        """
        Stoppt den Kratos Formoptimierungsprozess

        Args:
            None

        Returns:
            None
        """
        print("Stop Kratos optimization run.")
        print("#"*80)
        exit(0)

    def ReadValue(self,
                  identifier: str,
                  current_iteration_path: os.PathLike,
                  current_design: Kratos.Model(),
                  iteration: int) -> float:
        """
        Gibt den Wert der Antwortfunktion mit dem Namen des
        'identifier' aus, welcher in der .json der Kratos
        Optimierungsparameter spezifiziert wurde.

        Z.b. Masse als Zielfunktion in json:
        "objectives" : [{
            "identifier" : "mass",
            "type"       : "minimization",
            "analyzer"   : "optistruct"
        }]

        Args:
            identifier:                 einzigartiger Name der Antwortfunktion
            current_iteration_path:     Pfad des derzeitigen Iterationsschrittes
            current_design:             Derzeitiges Kratos Model
            iteration:                  Zahl des derzeitigen Iterationsschrittes (erste Iteration = 1)

        Returns:
            value:                      Wert der Antwortfunktion
        """
        NotImplementedError("ExternalSolver:: ReadValue has to be implemented!")

    def ReadShapeGradient(self,
                          identifier: str,
                          current_iteration_path: os.PathLike,
                          current_design: Kratos.Model(),
                          iteration: int) -> dict:
        """
        Gibt die Formsensitivitäten der Antwortfunktion mit dem Namen des
        'identifier' aus, welcher in der .json der Kratos
        Optimierungsparameter spezifiziert wurde.

        Sensitivitäten sind in folgendem Format zu übergeben:

        - shape_gradient:   python Dictionary
                            Schlüssel:  Knoten Id
                            Werte:      Formsensitivität in globalen kartesischen Koordinaten
                            {Knoten Id: [dF/dX, dF/dY, dF/dZ]}
                            {int: [float, float, float]}

        Args:
            identifier:                 einzigartiger Name der Antwortfunktion
            current_iteration_path:     Pfad des derzeitigen Iterationsschrittes
            current_design:             Derzeitiges Kratos Model
            iteration:                  Zahl des derzeitigen Iterationsschrittes (erste Iteration = 1)

        Returns:
            shape_gradient:             Formsensitivitäten als python Dict
        """
        NotImplementedError("ExternalSolver:: ReadShapeGradient has to be implemented!")

    def ReadThicknessGradient(self,
                              identifier: str,
                              current_iteration_path: os.PathLike,
                              current_design: Kratos.Model(),
                              iteration: int) -> dict:
        """
        Gibt die Schalendickensitivitäten der Antwortfunktion mit dem Namen des
        'identifier' aus, welcher in der .json der Kratos
        Optimierungsparameter spezifiziert wurde.

        Sensitivitäten sind in folgendem Format zu übergeben:

        - thickness_gradient:   python Dictionary
                                Schlüssel:  Knoten Id
                                Werte:      Schalendickensensitivität
                                {Knoten Id: dF/dT}
                                {int: float}

        Args:
            identifier:                 einzigartiger Name der Antwortfunktion
            current_iteration_path:     Pfad des derzeitigen Iterationsschrittes
            current_design:             Derzeitiges Kratos Model
            iteration:                  Zahl des derzeitigen Iterationsschrittes (erste Iteration = 1)

        Returns:
            thickness_gradient:         Schalendickensitivitäten als python Dict
        """
        NotImplementedError("ExternalSolver:: ReadThicknessGradient has to be implemented!")
