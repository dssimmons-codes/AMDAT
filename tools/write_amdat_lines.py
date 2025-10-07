#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on Unknown

this file should read in system.data/input_000.data
and write out the chemistry data required for AMDAT input
note system.data has excessive atom types

@author: forrest8
"""

# import os
from collections import Counter
import sys

mol_i = []  # molecule index
atom_type = []
# in automation, current working directory is a mess
# dir_path = os.path.dirname(os.path.realpath(__file__))
# print(dir_path) #magic
dir_input = "./lt_files/system.data"
# use cmd arg if given
if len(sys.argv) > 1:
  dir_input = sys.argv[1]
# dir_output = 'amdat_chemistry.inp'
with open(dir_input) as f:
    atom_switch = True
    atom_count = 10  # arbitrary int
    while atom_switch:
        split_line = f.readline().split()
        if len(split_line) > 1 and split_line[1] == "atoms":
            atom_count = int(split_line[0])
        if atom_switch and len(split_line) > 0 and split_line[0] == "Atoms":
            atom_switch = False  # atom header
    f.readline()  # empty
    for i in range(atom_count):
        m_t = f.readline().split()[1:3]
        mol_i.append(int(m_t[0]))
        atom_type.append(int(m_t[1]))
# done
molecule_count = mol_i[-1]
all_types = sorted(set(atom_type))

counter_init = Counter({x: 0 for x in all_types})  # inherits order/sortedness

poly_counters = [counter_init.copy() for i in range(molecule_count)]
print(poly_counters)
# mtc[0][int(mol_type[0][1])] += 1
# print(mtc[0])
for m_t in zip(mol_i, atom_type):
    poly_counters[m_t[0] - 1][m_t[1]] += 1

# remains to format for amdat

poly_line = " ".join(["poly" + str(i + 1) + " 1" for i in range(molecule_count)])

type_line = " ".join(map(str, all_types))

type_lines = [" ".join(map(str, poly.values())) for poly in poly_counters]

amdat_template = open("template_amdat.inp").read().splitlines()
template_head = amdat_template[0:4]
template_tail = amdat_template[8:]

with open("amdat.inp", "w") as f:
    f.write("\n".join(template_head) + "\n")
    # list(map(lambda x: x+"\n",amdat_template))
    # that's just worse isn't it
    f.write(poly_line + "\n")
    f.write(type_line + "\n")
    f.write("\n".join(type_lines) + "\n\n")
    f.write("\n".join(template_tail) + "\n")

"""
with open("amdat_lines.txt", "w") as f:
    f.write(poly_line + "\n")
    f.write(type_line + "\n")
    f.write("\n".join(type_lines) + "\n")
"""

