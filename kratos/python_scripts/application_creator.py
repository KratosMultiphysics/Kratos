import os
import re
import sys
import shutil


def ToUpperFromCamel(appCamel):
    ''' Converts a Camel-Case string into a upercase snake_case string '''

    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', appCamel)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).upper()


def ToLowerFromCamel(appCamel):
    ''' Converts a Camel-Case string into a lowercase snake_case string'''

    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', appCamel)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def CheckNameAvail(appsdir, appname):
    ''' Checks that the application does not exists

    Raises
    ------
    Error :
        If the application "appname" already exists

    '''

    appDir = appsdir + appname + "Application"

    if os.path.exists(appDir):
        raise Exception('Application "{}" already exsists'.format(appname))


def GetApplicationsDirectory():
    ''' Return the path to the applications directory '''

    return os.path.dirname(os.path.realpath(__file__)) + "/../../applications/"


def BuildApplication(apps, appCamel, appCaps, appLow):
    ''' Generates the application by copying it from the template app '''

    root = os.path.dirname(os.path.realpath(__file__)) + "/"

    tpldir = root + "../templates/" + "template_application"
    appdir = apps + appCamel + 'Application'

    # Create the directory
    shutil.copytree(tpldir, appdir, symlinks=False, ignore=None)

    # Copy all files from the template application
    filesToRename = {
        "template_application.cpp.in": appLow+"_application.cpp.in",
        "template_application.h.in": appLow+"_application.h.in",
        "TemplateApplication.py.in": appCamel+"Application.py.in",
        "tests/test_TemplateApplication.py.in":
            "tests/test_"+appCamel+"Application.py.in",
        "custom_python/test_python_application.cpp.in":
            "custom_python/"+appLow+"_python_application.cpp.in"
    }

    # Rename the files
    for key, value in filesToRename.items():
        os.rename(
            os.path.join(appdir, key),
            os.path.join(appdir, value)
        )

    # Replace the tokens in the app files
    for root, subfolder, files in os.walk(appdir):
        for f in files:

            src = os.path.join(root, f)
            fileName = src.split(".")
            dst = '.'.join(src.split(".")[0:len(fileName)-1])

            with open(src) as srcFile, open(dst, 'w+') as dstFile:
                for l in srcFile:

                    l = l.replace("@{APP_NAME_CAMEL}", appCamel)
                    l = l.replace("@{APP_NAME_CAPS}", appCaps)
                    l = l.replace("@{APP_NAME_LOW}", appLow)
                    l = l.replace("@{APP_DEPEND_LIST}", "")

                    dstFile.write(l)

            os.remove(src)


def main():

    # Fetch the applications directory
    appDir = GetApplicationsDirectory()

    # Read the application name and generate Camel, Caps and Low
    appCamel = sys.argv[1]
    appCaps = ToUpperFromCamel(appCamel)
    appLow = appCaps.lower()

    # Check that no application with this name currently exists
    CheckNameAvail(appDir, appCamel)

    # Generate the app
    BuildApplication(appDir, appCamel, appCaps, appLow)

    # Remind the user to change other files
    print("Your application has been generated in: applications/{}Application".format(appCamel))

if __name__ == "__main__":
    main()
