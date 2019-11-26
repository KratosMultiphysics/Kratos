from subprocess import Popen, PIPE

#----------------------------------------------------------------------------
# Get version string from git
#
# Author: Douglas Creager <dcreager@dcreager.net>
# http://dcreager.net/2010/02/10/setuptools-git-version-numbers/
#
# PEP 386 adaptation from
# https://gist.github.com/ilogue/2567778/f6661ea2c12c070851b2dfb4da8840a6641914bc
#----------------------------------------------------------------------------
def call_git_describe(abbrev=4):
    try:
        p = Popen(['git', 'describe', '--abbrev=%d' % abbrev],
                  stdout=PIPE, stderr=PIPE)
        p.stderr.close()
        line = p.stdout.readlines()[0]
        return line.strip().decode('utf8')

    except:
        return None


def pep386adapt(version):
    # adapt git-describe version to be in line with PEP 386
    parts = version.split('-')
    if len(parts) > 1:
        parts[-2] = 'post'+parts[-2]
        version = '.'.join(parts[:-1])
    return version


def git_version(abbrev=4):
    # First try to get the current version using "git describe".
    version = call_git_describe(abbrev)

    # If we still don't have anything, that's an error.

    if version is None:
        raise ValueError("Cannot find the version number!")

    #adapt to PEP 386 compatible versioning scheme
    version = pep386adapt(version)

    return version
