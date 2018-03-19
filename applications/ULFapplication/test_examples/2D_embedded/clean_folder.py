from glob import glob
import os

files_list = glob('*.time')
for f in files_list:
    os.remove(f)

files_list = glob('*.cond')
for f in files_list:
    os.remove(f)

files_list = glob('*.elem')
for f in files_list:
    os.remove(f)

files_list = glob('*.init')
for f in files_list:
    os.remove(f)

files_list = glob('*.node')
for f in files_list:
    os.remove(f)

files_list = glob('*.post.lst')
for f in files_list:
    os.remove(f)

files_list = glob('*.prop')
for f in files_list:
    os.remove(f)
