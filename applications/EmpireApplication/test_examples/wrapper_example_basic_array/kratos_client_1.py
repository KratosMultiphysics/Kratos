from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
# EMPIRE
import KratosMultiphysics.EmpireApplication as KratosEmpire
from ctypes import cdll
import os
import ctypes as ctp

print("This is kratos_client_1")

print("Starting to initialize Empire")
import empire_wrapper
print("Import Successfull")
empire = empire_wrapper.EmpireWrapper(echo_level=2)
print("Wrapper Created")
empire.Connect("kratos_client_1.xml")

array_to_send = [1.5,2.6,3,4,5]
empire.SendArray("array_1", array_to_send)
print("Sent Array:", array_to_send)

received_array = empire.ReceiveArray("array_2", 10)
print("Received Array:", received_array)

empire.Disconnect()