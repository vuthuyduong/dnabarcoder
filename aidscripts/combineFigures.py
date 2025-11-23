#!/usr/bin/env python
# -*- coding: utf-8 -*-
# FILE: combineFigures.py
# AUTHOR: Duong Vu
# CREATE DATE: 08 oct 2020

import os, argparse
import matplotlib.pyplot as plt
from PIL import Image

parser=argparse.ArgumentParser(prog='combineFigures.py', 
							   usage="%(prog)s [options] -i listof4figures -o outputprefix",
							   description='''Script that plot predictions from different predictionfiles for a taxonomic clade''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )
parser.add_argument('-i','--input', required=True, help='the list of images, separated by commas.')
parser.add_argument('-o','--out', required=True, help='The output filename.')

args=parser.parse_args()
imagefilenames=args.input
output=args.out

images=imagefilenames.split(",")
# --- Load images ---
imgA = Image.open(images[0])
imgB = Image.open(images[1])
imgC = Image.open(images[2])
imgD = Image.open(images[3])

# --- Create a 2x2 figure ---
fig, axes = plt.subplots(2, 2, figsize=(8,6))
plt.subplots_adjust(left=0, right=1, bottom=0, top=1,
                    wspace=0.01, hspace=0.01)

# Remove axes borders
for ax in axes.flatten():
    ax.axis("off")

# Place images
axes[0, 0].imshow(imgA)
axes[0, 1].imshow(imgB)
axes[1, 0].imshow(imgC)
axes[1, 1].imshow(imgD)

# Save the combined figure
#plt.tight_layout()
plt.savefig(output + ".pdf", dpi=600, bbox_inches="tight", pad_inches=0)
plt.savefig(output + ".png", dpi=600, bbox_inches="tight", pad_inches=0)
plt.close()
