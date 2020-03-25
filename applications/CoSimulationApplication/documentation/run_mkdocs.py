import glob
import os
import shutil

"""
README
------

- execute from this directory to update documentation
- required software:
    - mkdocs --> install with: pip install mkdocs
    - material theme --> install with: pip install mkdocs-material
"""

# clean docs folder and add coconuts
shutil.rmtree('docs')

os.mkdir('docs')
os.mkdir('docs/images')
os.mkdir('docs/assets')
os.mkdir('docs/assets/images')

shutil.copy('logo.png', 'docs/images/logo.png')
shutil.copy('favicon.ico', 'docs/assets/images/favicon.ico')

# find all MarkDown files in CoCoNuT
files = glob.glob('../**/*.md', recursive=True)

# check for duplicate filenames
filenames = []
for file in files:
    filenames.append(file.split('/')[-1])

for i, filename in enumerate(filenames):
    if filenames.count(filename) > 1:
        print(f'WARNING - duplicate file "{files[i]}"')

# copy all MarkDown files to docs folder
for file in files:
    shutil.copy(file, 'docs/')

# check if all files are mentioned in nav
unused = []
used = False
for filename in filenames:
    with open('mkdocs.yml', 'r') as file:
        for line in file:
            if filename in line:
                used = True
                break
    if not used:
        unused.append(filename)
    used = False
for file in unused:
    print(f'WARNING - file "{file}" is not used in mkdocs.yml')

# build static website
print('\n')
os.system('mkdocs build --clean')

# deploy website (need admin privileges on GitHub)
os.system('mkdocs gh-deploy')
