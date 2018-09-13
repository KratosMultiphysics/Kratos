#!/bin/bash

# This file is used to make a package of kratos in linux. It has been optimized to run with
# the 16.04 LTS of ubuntu

KRATOS_MAJOR="5" 																# Kratos major
KRATOS_MINOR="3" 																# Kratos minor
KRATOS_PATCH="0" 																# Kratos minor
KRATOS_ARCHY="64"
KRATOS_TDATE="XX-XX-XX"

KRATOS_RLS="$HOME/KratosRelease"
KRATOS_SRC="$HOME/KratosSource"
KRATOS_CMP="$HOME/KratosInstall"
KRATOS_AUX="$HOME/KratosLibsForPack64"

GIDINT_SRC="$HOME/GidInterface"
GIDINT_KTS=""

DEPLOY_DIR="$HOME/KratosDeploy"

# Obtain the Git hash for the current version
GIT_HASH=`git --git-dir="${KRATOS_SRC}/.git" rev-parse --short HEAD`
GIT_BRAN=`git --git-dir="${KRATOS_SRC}/.git" rev-parse --abbrev-ref HEAD`
GIT_NUMB="${GIT_HASH}(${GIT_BRAN})"

# Clean up old releases and files in the stage folder to prevent problems
rm -rf ${KRATOS_RLS}/*

# Copy Kratos, Kratos.gid and additional libs to the local release dir
cp -r ${GIDINT_SRC}/* ${KRATOS_RLS}
cp -r ${KRATOS_CMP}/* ${KRATOS_RLS}/kratos.gid

# Remove simlinks to system libs that are from not default packages
rm ${KRATOS_RLS}/kratos.gid/libs/libtrilinos*
rm ${KRATOS_RLS}/kratos.gid/libs/libmetis*

# Copy the extra needed libs
cp -r ${KRATOS_AUX}/* ${KRATOS_RLS}/kratos.gid

# Create a tar
tar -czf "${KRATOS_RLS}/kratos-${KRATOS_MAJOR}.${KRATOS_MINOR}.${KRATOS_PATCH}-${GIT_NUMB}-linux-${KRATOS_ARCHY}.tar.gz" -C "${KRATOS_RLS}" kratos.gid

cp "${KRATOS_RLS}/kratos-${KRATOS_MAJOR}.${KRATOS_MINOR}.${KRATOS_PATCH}-${GIT_NUMB}-linux-${KRATOS_ARCHY}.tar.gz" "${DEPLOY_DIR}/kratos-${KRATOS_MAJOR}.${KRATOS_MINOR}.${KRATOS_PATCH}-${GIT_NUMB}-linux-${KRATOS_ARCHY}.tar.gz"
