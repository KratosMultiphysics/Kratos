from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import os
#import kratos core and applications
import KratosMultiphysics

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

#Utility class to manage json settings
class JsonSettingsUtility(object):

    @staticmethod
    # delete not matching settings supplied settings and print a warning
    def CheckAndFixNotMatchingSettings(supplied_settings, expected_settings):
        """Check and fix not matching settings from origin to destination.
        """
        update_needed = False
        # check if the parameter exists in the expected results
        for name, supplied_value in supplied_settings.items():
            if not expected_settings.Has(name): #  transfer value.
                print(" INPUT PARAMETER NOT CONSIDERED: ["+name+"] " )
                supplied_settings.RemoveValue(name)
                update_needed = True

        if( update_needed ):
            print(" ::[PLEASE consider to UPDATE the INPUT FILE to run the case with the GIVEN SETTINGS]:: ")

        # check if the parameter has the same type as expected
        for name, expected_value in expected_settings.items():
            if supplied_settings.Has(name): #  transfer value.
                supplied_value = supplied_settings[name]
                if JsonSettingsUtility.CheckMatchingValueType(supplied_value,expected_value):
                    if expected_value.IsArray() and supplied_value.IsArray():
                        if expected_value.size() == supplied_value.size():
                            for i in range(expected_value.size()):
                                if JsonSettingsUtility.CheckMatchingValueType(supplied_value[i],expected_value[i]):
                                    if expected_value[i].IsSubParameter() and supplied_value[i].IsSubParameter():
                                        JsonSettingsUtility.CheckAndFixNotMatchingSettings(supplied_value[i], expected_value[i])
                                else:
                                    raise Exception('Unsupported sub parameter type: ' + name)
                    elif expected_value.IsSubParameter() and supplied_value.IsSubParameter():
                        JsonSettingsUtility.CheckAndFixNotMatchingSettings(supplied_value, expected_value)
                else:
                    raise Exception('Unsupported parameter type: ' + name)


    @staticmethod
    # transfer from origin to destination: custom settings in derived classes transferred and deleted from the supplied settings.
    def TransferMatchingSettingsToDestination(origin_settings, destination_settings):
        """Transfer matching settings from origin to destination.

        If there is any name/value in the origin settings matching with the destination settings,
        then the setting value is assigned to the destination, and deleted from the origin.

        """
        #print("start",origin_settings.PrettyPrintJsonString())
        for name, destination_value in destination_settings.items():
            if origin_settings.Has(name): # Validate and transfer value.
                origin_value = origin_settings[name]
                if JsonSettingsUtility.CheckAndTransferMatchingValueType(origin_value,destination_value):
                    if destination_value.IsArray() and origin_value.IsArray():
                        if destination_value.size() != origin_value.size():
                            raise Exception('len("' + name + '") != ' + str(destination_value.size()))
                        for i in range(destination_value.size()):
                            if JsonSettingsUtility.CheckAndTransferMatchingValueType(origin_value[i],destination_value[i]):
                                if destination_value[i].IsSubParameter() and origin_value[i].IsSubParameter():
                                    JsonSettingsUtility.TransferMatchingSettingsToDestination(origin_value[i], destination_value[i])
                                    if len(origin_value[i].items()) != 0:
                                        raise Exception('Json settings not found in default settings: ' + origin_value[i].PrettyPrintJsonString())
                            else:
                                raise Exception('Unsupported sub parameter type: ' + name)
                    elif destination_value.IsSubParameter() and origin_value.IsSubParameter():
                        JsonSettingsUtility.TransferMatchingSettingsToDestination(origin_value, destination_value)
                        if len(origin_value.items()) != 0:
                            raise Exception('Json settings not found in default settings: ' + origin_value.PrettyPrintJsonString())
                else:
                    raise Exception('Unsupported parameter type: ' + name)
                origin_settings.RemoveValue(name)
        #print("end",origin_settings.PrettyPrintJsonString())
        #print("result",destination_settings.PrettyPrintJsonString())

    @staticmethod
    # checks and transfers mathing values (except array and subparameters)
    def CheckAndTransferMatchingValueType(origin_value, destination_value):
        if destination_value.IsDouble() and origin_value.IsDouble():
            destination_value.SetDouble(origin_value.GetDouble())
            return True
        elif destination_value.IsInt() and origin_value.IsInt():
            destination_value.SetInt(origin_value.GetInt())
            return True
        elif destination_value.IsBool() and origin_value.IsBool():
            destination_value.SetBool(origin_value.GetBool())
            return True
        elif destination_value.IsString() and origin_value.IsString():
            destination_value.SetString(origin_value.GetString())
            return True
        elif destination_value.IsArray() and origin_value.IsArray():
            return True
        elif destination_value.IsSubParameter() and origin_value.IsSubParameter():
            return True
        else:
            return False

    @staticmethod
    # checks matching values
    def CheckMatchingValueType(origin_value, destination_value):
        if destination_value.IsDouble() and origin_value.IsDouble():
            return True
        elif destination_value.IsInt() and origin_value.IsInt():
            return True
        elif destination_value.IsBool() and origin_value.IsBool():
            return True
        elif destination_value.IsString() and origin_value.IsString():
            return True
        elif destination_value.IsArray() and origin_value.IsArray():
            return True
        elif destination_value.IsSubParameter() and origin_value.IsSubParameter():
            return True
        else:
            return False
