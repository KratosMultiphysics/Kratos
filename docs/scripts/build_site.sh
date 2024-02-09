#!/bin/bash

if [ ! -d ../_site/Page_files ]; then
    if [ -d ../_site ]; then
        rm -rf ../_site
    fi

    git clone --depth 1 --filter=blob:none --sparse https://github.com/KratosMultiphysics/Documentation.git ../_site
    cd ../_site
    git sparse-checkout init --cone
    git sparse-checkout set Page_files

    rm ../_site/*
fi

cp -r ../pages ../_site/Page_files
cp ../_config.yml ../_site/Page_files
if [ ! -d "../_site/Page_files/_data" ] || [ "$1" == "build_menus" ]; then
    cd ..
    python3 process_pages.py -t local
    cd scripts
fi

if [ -d "../_data" ]; then
    if [ -d ../_site/Page_files/_data ]; then
        rm -r ../_site/Page_files/_data
    fi
    mv ../_data ../_site/Page_files/_data
fi

cd ../_site/Page_files
jekyll build