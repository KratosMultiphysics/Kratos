#!/bin/bash

wget --convert-links -r http://127.0.0.1:4000/
rm -r packaged_site
mv 127.0.0.1:4000 packaged_site