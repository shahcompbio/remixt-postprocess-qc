import argparse

import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import remixt
import remixt.cn_plot as cn_plot

def parse_args():
    description = "Plot ReMixT-Postprocessing QC plots"
    p = argparse.ArgumentParser(description=description)
    p.add_argument('-i', '--input', help="remixt_cn path", required=True)
    p.add_argument('-o', '--output', help="output QC png path", required=True)
    args = p.parse_args()
    return args

def add_legend(fig, label_colors, marker='_'):
    ax = fig.get_axes()[0]
    handle = [plt.plot([], [],
              color=label_colors[label], marker=marker, ms=4, ls="")[0] 
              for label in label_colors]
    legend = ax.legend(handles=handle, labels=label_colors.keys(), title="labels")
    legend.get_frame().set_alpha(0.4)
    ax.add_artist(legend);

if __name__ == "__main__":
    args = parse_args()
    data = pd.read_csv(args.input, sep='\t', dtype={'chromosome': 'str'})
    data['chromosome'] = data['chromosome'].str.replace('chr', '')
    p = PdfPages(args.output)

    # fig 1
    fig = plt.figure(figsize=(8,8))
    box = matplotlib.transforms.Bbox([[0., 0.], [1., 1.]])
    transform = matplotlib.transforms.BboxTransformTo(box)
    _ = cn_plot.plot_cnv_scatter_density(fig, transform, data, major_col='major_raw', minor_col='minor_raw')
    ax1, _, ax3, ax4 = fig.get_axes()
    legend = ax1.get_legend()
    legend.set_bbox_to_anchor((0.8, 0.4))
    p.savefig(fig)

    # fig 2
    fig = plt.figure(figsize=(15, 3))
    box = matplotlib.transforms.Bbox([[0., 0.], [1., 1.]])
    transform = matplotlib.transforms.BboxTransformTo(box)
    chromosomes = [str(a) for a in list(range(1, 23)) + ['X', 'Y']]
    _ = cn_plot.plot_cnv_genome_density(fig, transform, data, chromosomes=chromosomes, scatter=True)
    label_colors = {'major_raw': 'tab:red', 'minor_raw':'tab:blue'}
    add_legend(fig, label_colors, marker='o')
    p.savefig(fig)

    # fig 3
    fig, ax = plt.subplots(figsize=(11.8, 2))
    _ = cn_plot.plot_cnv_genome(
        plt.gca(), data, major_col='major_raw_e', minor_col='minor_raw_e',
        chromosomes=chromosomes, scatter=False, do_fill=True, maxcopies=8)
    label_colors = {'major_raw_e': 'tab:red', 'minor_raw_e':'tab:blue'}
    add_legend(fig, label_colors)
    p.savefig(fig)

    # fig 4
    fig, ax = plt.subplots(figsize=(11.8, 2))
    _ = cn_plot.plot_cnv_genome(
        plt.gca(), data, major_col='major_1', minor_col='minor_1',
        chromosomes=chromosomes, scatter=False, do_fill=True, maxcopies=8)
    label_colors = {'major_1': 'tab:red', 'minor_1':'tab:blue'}
    add_legend(fig, label_colors)
    p.savefig(fig)

    # close
    p.close()
