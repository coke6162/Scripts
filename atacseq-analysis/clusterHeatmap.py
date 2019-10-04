#!/usr/bin/python

## Script for making heatmap from giggle score matrix
## Date: 10 May 2019

## Example usage:
## cd /Users/coke6162/repos/scripts/atacseq-analysis
## /opt/python/3.6.3/bin/python3.6 \ 
## ./clustermap.py \
## -i /Users/coke6162/project20xx_atacseq/giggle/allRep/giggleMatrix_allRep.tab \
## -x 15 -y 40 -m ward --no_col_cluster \
## -o /Users/coke6162/project20xx_atacseq/giggle/allRep/giggleMatrix_allRep.png

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import math
from optparse import OptionParser
from scipy.spatial.distance import pdist,squareform

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

parser = OptionParser()

parser.add_option("-i",
                  "--input",
                  dest="input_file",
                  help="GIGGLE results matrix of combo scores")

parser.add_option("-o",
                  "--output",
                  dest="output_file",
                  help="Output file name")

parser.add_option("-x",
                  "--x_size",
                  dest="x_size",
                  type="int",
                  default=10,
                  help="Figure x size (Default 10)")

parser.add_option("-y",
                  "--y_size",
                  dest="y_size",
                  type="int",
                  default=30,
                  help="Figure y size (Default 30)")

parser.add_option("--col_cluster", 
                  action="store_true", 
                  dest="col_cluster",
                  help="Cluster columns")
parser.add_option("--no_col_cluster", 
                  action="store_false", 
                  dest="col_cluster",
                  help="Don't cluster columns")

parser.add_option("-m",
                  "--method",
                  dest="cluster_method",
                  type="string",
                  default="ward",
                  help="Cluster method (Default 'ward'). Other options: single, complete, average, weighted, centroid, median")

(options, args) = parser.parse_args()

if not options.input_file:
    parser.error('Input file not given')
if not options.output_file:
    parser.error('Output file not given')


data = pd.read_table(options.input_file, sep='\t', index_col=0)
# print(data) 

# Make a heatmap with hierarchial clustering and save output figure
import seaborn as sns; sns.set(color_codes=True)

# choose colour scheme
# copying Ryan's colour scheme from giggle paper
from matplotlib import colors as mcolors
from matplotlib.colors import Normalize

_seismic_data = ( (0.0, 0.0, 0.3), 
                  (0.0, 0.0, 1.0),

                  (1.0, 1.0, 1.0),

                  (1.0, 0.0, 0.0),
                  (0.5, 0.0, 0.0))

hm = mcolors.LinearSegmentedColormap.from_list( \
        name='red_white_blue', \
        colors=_seismic_data, N=256)

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

fig = plt.figure()
clustermap = sns.clustermap(data, cmap=hm, norm = MidpointNormalize(midpoint=0), robust=False, metric='euclidean', method=options.cluster_method, figsize=(options.x_size, options.y_size), col_cluster=options.col_cluster)
plt.savefig(options.output_file,bbox_inches='tight')