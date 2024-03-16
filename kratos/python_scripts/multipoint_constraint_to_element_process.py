import KratosMultiphysics

MultipointConstraintToElementProcess = KratosMultiphysics.MultipointConstraintToElementProcess

def Factory(parameters: KratosMultiphysics.Parameters,
            model: KratosMultiphysics.Model) -> MultipointConstraintToElementProcess:
    return MultipointConstraintToElementProcess(model, parameters["Parameters"])
