# SIA-IGA Optimierungs-Workflow (Reminder)  ####CHECKLEO

## Ziel
Verwendung der SystemIdentificationApplication (SIA) zur Optimierung von Aktuierungen (z.B. α) über Displacement-Sensoren mit adjungiertem Finite-Differenzen-Element aus der IGA Application.

## Schritte

1. **Input-Files kopieren**
   - Kopiere alle relevanten Input-Files aus:
     `kratos/applications/SystemIdentificationApplication/tests/auxiliary_files/system_identification/`
     in deinen Arbeitsordner (z.B. `092_patch1x1_sup-line-x_GD_Sub1-Opt`).
   - Typische Files: `sensor_data.json`, `primal_material_properties.json`, `adjoint_project_parameters.json`, ggf. weitere Material-/Parameterdateien.

2. **Input-Files anpassen**
   - Passe die Sensoren, Materialdaten und Pfade in den JSON-Files an deine Geometrie und Optimierungsaufgabe an.
   - Displacement-Sensoren reichen für die Steuerung der Aktuierungen.

3. **Adjoint-Element**
   - Das adjoint Element für Finite-Differenzen ist in der IGA Application implementiert.
   - Die Sensitivitätsberechnung (Gradient) erfolgt per Finite Differenzen im Element (ggf. anpassen).

4. **SIA-Analyse starten**
   - Starte die SIA mit einem Python-Skript, z.B.:
     ```python
     import KratosMultiphysics
     import KratosMultiphysics.SystemIdentificationApplication as SIA
     from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers import sensor_sensitivity_static_analysis
     import json
     with open("adjoint_project_parameters.json") as f:
         parameters = KratosMultiphysics.Parameters(f.read())
     model = KratosMultiphysics.Model()
     analysis = sensor_sensitivity_static_analysis(model, parameters)
     analysis.Run()
     ```

5. **Optimierungsschleife**
   - Baue eine Schleife um die SIA-Analyse:
     - Führe SIA aus → erhalte Sensitivitäten/Gradienten → update Aktuierung (α) → schreibe neue Input-Files → wiederhole.
   - Kein Custom-Sensor nötig, solange nur Displacements kontrolliert werden.

## Hinweise
- Die eigentliche Optimierungslogik (Update der Aktuierung) liegt außerhalb der SIA und wird in einem eigenen Skript umgesetzt.
- Die SIA liefert die Sensitivitäten/Gradienten, die für das Update benötigt werden.
- Die Finite-Differenzen-Logik im adjoint Element kann bei Bedarf angepasst werden.

---
Letzte Aktualisierung: 2025-10-02
