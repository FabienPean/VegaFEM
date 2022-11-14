from itertools import count
import os
import re
import csv
from pathlib import Path


def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return root.split(os.path.sep)[-1]

rootdir = 'libraries'

directory_list = list()
for subdir, dirs, files in os.walk(rootdir):
    for name in dirs:
        directory_list.append(name)

subfolders = [f for f in Path(rootdir).iterdir() if f.is_dir()]

public_dependencies = dict()
for subdir, dirs, files in os.walk(rootdir):
    dirname = subdir.split(os.path.sep)[-1]
    dependencies = list()
    for file in files:
        if file.endswith(".h") or file.endswith(".hpp"):
            file_path = os.path.join(subdir, file)
            #print(file_path)
            with open(file_path) as f:
                txt = f.read()
                txt = re.sub("#ifndef[ \t]*?(?!VEGAFEM_)_*(\w*?)_*?[ \t]*?\n(#define[ \t]*?\w*?)[ \t]*?\n",'#ifndef VEGAFEM_\\1\n\\2\n', txt, flags=re.DOTALL, count=1)
                # matches: "#ifndef"  there, then capture any word (ascii characters and underscore) as many times as necessar
                # followed by either tab or space 0 or N times non-eager: [ \t]*? 
                # discards the regex result if it contains "VEGAFEM_" there: (?!VEGAFEM_)
                # matches all underscores: _*
                # captures all characters (ascii chars and underscore) 0 or N times non-eager: (\w*?)
                # matches all underscores followed by either tab or space 0 or N times non-eager until the first newline character encountered: _*?[ \t]*?\n
                # replaces the whole string match by '#ifndef VEGAFEM_\\1\n\\2\n' wherein \\1 and \\2 are, respectively, the first and second capture groups
                # then captures the second line according to: (#define[ \t]*?\w*?)[ \t]*?\n
                # re.DOTALL makes * catches absolutely every symbol (including new line and others)
                # count=1 for one occurrence
                txt = re.sub("(#ifndef[ \t]*?\w*?)[ \t]*?\n#define[ \t]*?(?!VEGAFEM_)_*(\w*?)_*?[ \t]*?\n",'\\1\n#define VEGAFEM_\\2\n', txt, flags=re.DOTALL, count=1)
            with open(file_path, "w") as f:
                f.write(txt)

