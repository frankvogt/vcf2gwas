import os
import shutil

"""
Copyright (C) 2021, Frank Vogt

This file is part of vcf2gwas.

vcf2gwas is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vcf2gwas is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with vcf2gwas.  If not, see <https://www.gnu.org/licenses/>.
"""

def main():

    print("Copying example input files to current working directory..")
    try:
        os.mkdir("input")
    except:
        pass
    path = os.path.join(os.path.dirname(__file__), 'input')
    for file in os.listdir(path):
        shutil.copy(os.path.join(path, file), os.path.join('input', file))
    print("Copying README and LICENSE files..")
    shutil.copy(os.path.join(os.path.dirname(__file__), 'README.pdf'), 'README.pdf')
    shutil.copy(os.path.join(os.path.dirname(__file__), 'LICENSE'), 'LICENSE')
