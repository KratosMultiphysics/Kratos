# Importing the Kratos Library
import KratosMultiphysics as KM


def GetBDFIntegrationOrder(scheme_type):
    if scheme_type == "backward_euler":
        order = 1
    else:
        if scheme_type == "bdf":
            raise Exception('Wrong input for scheme type: "bdf"! Please append the order to the bdf-scheme, e.g. "bdf2"')
        # BDF schemes can be from 1 to 5 order, so in order to detect the integration order from the scheme_type we remove the "bdf" string, that is, if the user tells bdf3 only 3 will remain when we remove bdf which corresponds to the method of choice
        order = int(scheme_type.replace("bdf", ""))

    # Warning
    if (order > 2):
        KM.Logger.PrintWarning("BDF", "Order {}; constant time step must be considered".format(order))

    return order
