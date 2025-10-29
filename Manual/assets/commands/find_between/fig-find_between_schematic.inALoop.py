#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on Fri March 3 05:48:00 2023

@author: pierrekawak
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('tkagg')
import numpy as np
import sys

# figure settings
test             = 0
fontsize         = 8
labelsize        = 6
figsize          = (1.75, 1.75)
pt_to_in         = 72
neat_marker      = "^"
filled_marker    = "o"
filler_marker    = "s"
polymer_marker   = "D"
neat_color       = "r"
neat_color       = "#0B7A75"
neat_color       = "#5ADBFF"
neat_color       = "#70A9A1"
filled_color     = "#808080"
filled_color     = "#800080"
filler_color     = "b"
polymer_color    = "r"
other_color      = "c"

np.random.seed(15023)

add_to_rand = 1.5
center_1    = (0.0, 0.0)
center_2    = (3.0, 0.0)
numb_pos    = 2000000
#numb_pos    = 10000

ex_pos = (1.3, 0.88)

def get_diff(pos1, pos2):
  return tuple(map(lambda i, j: i - j, pos1, pos2))

def get_dist(pos1, pos2):
  return sum(tuple(map(lambda i, j: (i - j)**2, pos1, pos2)))

def get_angle(vec1, vec2):
  vec1_mag = sum(lambda i: i**2, vec1)
  vec2_mag = sum(lambda i: i**2, vec2)
  return np.arccos(sum(lambda i, j: i * j, vec1, vec2)/vec1_mag/vec2_mag)

def get_cos_angle(vec1, vec2):
#  vec1_mag = sum(tuple(map(lambda i: i**2, vec1)))
#  vec2_mag = sum(tuple(map(lambda i: i**2, vec2)))
  vec1_mag = np.linalg.norm(vec1)
  vec2_mag = np.linalg.norm(vec2)
  #return sum(tuple(map(lambda i, j: i * j, vec1/vec1_mag, vec2/vec2_mag)))
  return np.dot(vec1, vec2)/vec1_mag/vec2_mag

def check_condition_helper(hdist_1, hdist_2, cos_theta, dcut, tcut):
  if hdist_1 > dcut:
    return 0
  if hdist_2 > dcut:
    return 0
  if hdist_1 + hdist_2 > dcut:
    return 0
  if cos_theta > tcut:
    return 0
  return 1

def check_condition(pos, dcut, tcut):
  hdist_1    = get_dist(pos, center_1)/2
  hdist_2    = get_dist(pos, center_2)/2
  vec1       = tuple(map(lambda i, j: i - j, center_1, pos))
  vec2       = tuple(map(lambda i, j: i - j, center_2, pos))
  cos_theta  = get_cos_angle(vec1, vec2)
#  print(pos, center_1, center_2)
#  print(vec1, vec2, cos_theta, np.degrees(np.arccos(cos_theta)))
#  sys.exit()
  return(check_condition_helper(hdist_1, hdist_2, cos_theta, dcut, tcut))

max_pos_x = max(center_1[0], center_2[0])
min_pos_x = min(center_1[0], center_2[0])
max_pos_y = max(center_1[1], center_2[1])+add_to_rand
min_pos_y = min(center_1[1], center_2[1])-add_to_rand

# get numb_pos random floats between min_pos and max_pos
init_pos_x = np.random.uniform(low=min_pos_x, high=max_pos_x, size=(numb_pos,))
init_pos_y = np.random.uniform(low=min_pos_y, high=max_pos_y, size=(numb_pos,))
positions = [ (x, y) for x, y in zip(init_pos_x, init_pos_y) ]

dist_cutoffs      = [ 4.0, 4.0 , 5.0, 5.0 ]
cos_theta_cutoffs = [ 0.0, -0.5, 0.0, 0.5 ]
is_plot_examples  = [ 1  , 0   , 0  , 0   ]

for dist_cutoff, cos_theta_cutoff, is_plot_example in zip(dist_cutoffs, cos_theta_cutoffs, is_plot_examples):
  # check_condition returns 0 if pos does not satisfy the particles_between algorithm and 1 otherwise
  is_between = [ check_condition(pos, dist_cutoff, cos_theta_cutoff) for pos in positions ]
  # filter init_pos of all is_between == 0
  pos_x = np.array([ x for x, i in zip(init_pos_x, is_between) if i ])
  pos_y = np.array([ y for y, i in zip(init_pos_y, is_between) if i ])

  fig, ax = plt.subplots(1, figsize=figsize, sharey=True, constrained_layout=test)
  ax.scatter(pos_x, pos_y, s=10, color=polymer_color, alpha=1.0)#0.1)
  ax.scatter([center_1[0], center_2[0]], [center_1[1], center_2[1]], s=1000, color=filler_color)#, alpha=alpha, label=label)

  # example position
  if is_plot_example:
    vec1   = get_diff(center_1, ex_pos)
    vec2   = get_diff(center_2, ex_pos)
    ax.scatter(ex_pos[0], ex_pos[1], s=10, color='k')
    ax.arrow(ex_pos[0], ex_pos[1], vec1[0], vec1[1], head_width=0.1, length_includes_head=1, color='k')
    ax.arrow(ex_pos[0], ex_pos[1], vec2[0], vec2[1], head_width=0.1, length_includes_head=1, color='k')
    from AngleAnnotation import AngleAnnotation
    am1 = AngleAnnotation(ex_pos, center_1, center_2, ax=ax, size=300, text=r"$\theta$", textposition = "inside", text_kw=dict(fontsize=10, xytext = (-5,73)))
    ax.text((ex_pos[0]+center_1[0])/2*0.70, (ex_pos[1]+center_1[1])/2*1.30, r"$v_{1}$", multialignment="center")
    ax.text((ex_pos[0]+center_2[0])/2*0.92, (ex_pos[1]+center_1[1])/2*1.30, r"$v_{2}$", multialignment="center")

  ax.set_xlim(center_1[0]-0.1, center_2[0]+0.1)
  ylims = (ax.get_xlim()[1]-ax.get_xlim()[0])/2
  ax.set_ylim(-ylims, ylims)
  #print("xlim:", ax.get_xlim())
  #print("ylim:", ax.get_ylim())
  print("dist_cutoff:", dist_cutoff, "cos_theta_cutoff:", cos_theta_cutoff, "theta_cutoff [deg]:", np.degrees(np.arccos(cos_theta_cutoff)), "remaining%", len(pos_x)/len(init_pos_x)*100)

  is_plot_legend = 0
  if is_plot_legend:
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color=polymer_color, lw=4),
                  Line2D([0], [0], color=filler_color, marker="o")]
    custom_names = ["Condition Satisfied", "Filler Center Bead"]

    ax.legend(custom_lines, custom_names, loc="lower right", fontsize=fontsize)

  ax.tick_params(axis='both', labelsize=labelsize, pad=1)
  plt.axis('off')
  if test == 0:
    plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0)

  out_file_name = os.path.basename(sys.argv[0]).split('.')[0]+"_dc"+str(dist_cutoff)+"_tc"+str(cos_theta_cutoff)+"_eg"+str(is_plot_example)+".png"
  fig.savefig(out_file_name, dpi=900)
  print("Saving to", out_file_name)
