---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3.7.12 64-bit (conda)
    language: python
    name: python3712jvsc74a57bd0dc0fe05456373cce17991fe9e2e9264df4f4d1e972d77814a9700f21c9e7a8e2
---

# SNV CN


## modules

```python
%load_ext autoreload
%autoreload 2
```

```python
import glob
```

```python
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import remixt.cn_plot as cn_plot
import wgs_analysis.refgenome as refgenome
refgenome.set_genome_version('hg38')

from wgs_qc_utils.utils.empty import empty_plot
import wgs_qc_utils.plotter.input_checker as input_checker
```

```python
import os
```

```python
import yaml
```

```python
from wgs_qc_utils.reader import read_remixt
from wgs_qc_utils.reader import read_variant_calls
from wgs_qc_utils.reader import parse_snv_cn
```

```python
def plot_scatter(pos, frac_cn, axis, ylim=None, logistic_y=False):
    if not isinstance(pos, pd.Series) and not isinstance(frac_cn, pd.Series):
        return empty_plot(axis, "Snv Cn", snv_cn=True)

    input_checker.check_input_is_valid([pos, frac_cn],
                                            [input_checker.CheckerTypes.NUMERIC,
                                             input_checker.CheckerTypes.FLOAT])
    # axis2 = axis.twinx()
    axis2 = axis
    if ylim != None: axis2.set_ylim(ylim)
    # axis2.set_yticklabels([])
    # axis2.set_yticks([])
    if logistic_y:
        squash_coeff = 0.15
        squash_f = lambda a: np.tanh(squash_coeff * a)
        frac_cn = squash_f(frac_cn)
    # else:
        # frac_cn_ul = int(np.quantile(frac_cn, 0.95)) + 1
        # yticklabels = np.arange(0, frac_cn_ul+0.1, 0.5)
        # axis2.set_ylim(0, frac_cn_ul)
        # axis2.set_yticks(yticklabels)
        # axis2.set_yticklabels(yticklabels)
        frac_cn = frac_cn[frac_cn >= 0]

    axis2.scatter(pos, frac_cn, color="black", s=3, marker="o", alpha=0.2)
    return axis2

def plot_hist(frac_cn, axis, ylim=None, logistic_y=False):
    if not isinstance(frac_cn, pd.Series):
        return empty_plot(axis, "Snv cn", snv_cn=True)
    if frac_cn.empty:
        axis.set_ylim(0, 8)
        axis.set_ylabel("SNV density")
        return axis

    input_checker.check_input_is_valid([frac_cn],
                                            [input_checker.CheckerTypes.FLOAT])
    if logistic_y:
        squash_coeff = 0.15
        squash_f = lambda a: np.tanh(squash_coeff * a)
        frac_cn = squash_f(frac_cn)
        yticks = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 20])
        yticks_squashed = squash_f(yticks)
        ytick_labels = [str(a) for a in yticks]
        axis.set_yticks(yticks_squashed)
        axis.set_yticklabels(ytick_labels)
        axis.set_ylim((-0.01, 1.01))
        axis.spines['left'].set_bounds(0, 1)
    else:
        # frac_cn_ul = int(np.quantile(frac_cn, 0.95)) + 1
        # yticklabels = np.arange(0, frac_cn_ul+0.1, 0.5)
        # axis.set_ylim(0, frac_cn_ul)
        # axis.set_yticks(yticklabels)
        # axis.set_yticklabels(yticklabels)
        # frac_cn = frac_cn[frac_cn <= frac_cn_ul]
        frac_cn = frac_cn[frac_cn >= 0]
    if ylim != None: axis.set_ylim(ylim)
    axis.hist(frac_cn, color="black", orientation="horizontal", bins=60)
    axis.set_ylabel("SNV density")
    return axis


def add_legend(fig, label_colors, marker='_'):
    ax = fig.get_axes()[0]
    handle = [plt.plot([], [],
              color=label_colors[label], marker=marker, ms=4, ls="")[0]
              for label in label_colors]
    legend = ax.legend(handles=handle, labels=label_colors.keys(), title="labels")
    legend.get_frame().set_alpha(0.4)
    ax.add_artist(legend);
```

```python
def get_nygc_snvs(sample):
    vcf_cols = ['chrom', 'pos', 'ID', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'normal', 'tumor']
    info_cols = ['chrom', 'pos', 'ref', 'alt', 'DP', 'VAF']
    nygc_snv_dir = '/juno/work/shah/users/chois7/tickets/hcmi/downloads/snvs'
    vcf_paths = glob.glob(f'{nygc_snv_dir}/{sample}*--*.vcf.gz')
    if len(vcf_paths) == 0:
        return get_wusl_snvs(sample)
    vcf_path = vcf_paths[0]
    vdf = pd.read_table(vcf_path, comment='#', names=vcf_cols)
    vdf['chrom'] = vdf['chrom'].str.replace('chr', '')
    tdf = vdf['tumor'].str.split(':', expand=True)
    tdf.columns = ['AD', 'DP', 'VAF']
    tdf[['DP', 'VAF']] = tdf[['DP', 'VAF']].astype(float)
    vdf['DP'] = tdf['DP']
    vdf['VAF'] = tdf['VAF']
    return vdf[info_cols]
```

```python
def get_wusl_snvs(sample):
    src_cols = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2', 't_depth', 't_alt_count']
    dst_cols = ['chrom', 'pos', 'ref', 'alt', 't_depth', 't_alt_count']
    snv_dir = '/juno/work/shah/users/chois7/tickets/hcmi/downloads/snvs'
    maf_path = f'{snv_dir}/{sample}.result.maf'
    vdf = pd.read_table(maf_path, comment='#')
    vdf.rename(columns=dict(zip(src_cols, dst_cols)), inplace=True)
    vdf['chrom'] = vdf['chrom'].str.replace('chr', '')
    return vdf[dst_cols]
```

# Data

```python
meta_path = '/juno/work/shah/users/chois7/tickets/remixt-pp-plot/resources/paths.WGS-QC.tsv'
meta = pd.read_table(meta_path)
# samples = meta['isabl_sample_id'].unique()
```

## plot

```python
# samples = [
#     'HCM-BROD-0210-C71-85R'
# ]

samples = [
    "HCM-BROD-0038-C41-86A",
    "HCM-BROD-0195-C71-85A",
    "HCM-BROD-0199-C71-85A",
    "HCM-BROD-0210-C71-85R",
    "HCM-BROD-0214-C71-85A",
    "HCM-BROD-0220-C15-85B",
    "HCM-BROD-0344-C25-85A",
    "HCM-BROD-0448-C15-06A",
    "HCM-BROD-0577-C15-06A",
    "HCM-BROD-0578-C15-06A",
    "HCM-BROD-0695-C71-02B",
    "HCM-CSHL-0058-C34-85C",
    "HCM-CSHL-0058-C34-86A",
    "HCM-CSHL-0081-C25-01A",
    "HCM-CSHL-0176-C25-01A",
    "HCM-CSHL-0374-C25-01A",
    "HCM-CSHL-0379-C18-85A",
    "HCM-CSHL-0459-C17-01B",
    "HCM-CSHL-0729-C18-85A",
    "HCM-CSHL-0807-C54-01A",
    "HCM-SANG-0284-C18-01A",
    "HCM-SANG-0289-C15-85A",
    "HCM-SANG-0314-C15-85A",
]
```

```python
ploidy = 2.69
purity = 0.98
```

```python
snv_cn['corrected_VAF'] = snv_cn['VAF_tumor'] / (snv_cn['total_raw_e'] / ploidy)
snv_cn['estimated_ccf'] = snv_cn['corrected_VAF'] / purity
```

```python
tmp = snv_cn.copy()
```

```python
tmp.columns
```

```python
tmp['
```

```python

```

```python
for sample in samples[1:2]:
    paths = meta[meta['isabl_sample_id']==sample]
    paths = paths[paths['result_type']=='command_log']
    assert paths.shape[0] == 1, paths
    path = paths['result_filepath'].iloc[0]
    out_dir, _ = os.path.split(path)
    inputs_path = f'{out_dir}/inputs.yaml'
    assert os.path.exists(inputs_path), inputs_path
    inputs = yaml.load(open(inputs_path))
    
    remixt_path = inputs['remixt']
    remixt = read_remixt.read(remixt_path)
    remixt['chrom'] = remixt['chrom'].str.replace('chr', '')
    remixt = remixt.dropna()
    
    somatic_calls_path = inputs['somatic_calls']
    somatic_calls = read_variant_calls.read(somatic_calls_path)
    somatic_calls['chrom'] = somatic_calls['chrom'].str.replace('chr', '')
    
    nsnvs = get_nygc_snvs(sample)
    somatic_calls = somatic_calls[somatic_calls.pos.isin(nsnvs.pos)]
    
    snv_copynumber = parse_snv_cn.parse(somatic_calls, remixt)
    
    cnv = remixt.copy()
    cnv.rename(columns={'chrom':'chromosome'}, inplace=True)
    cnv['chromosome_start'] = cnv.chromosome.str.upper().map(refgenome.info.chromosome_start)
    cnv['coord'] = cnv['start'] + cnv['chromosome_start']
    snv_cn = snv_copynumber.copy().dropna()
    snv_cn['chromosome_start'] = snv_cn.chrom.str.upper().map(refgenome.info.chromosome_start)
    snv_cn['coord'] = snv_cn['pos'] + snv_cn['chromosome_start']
    snv_cn['snv_cn'] = (snv_cn['VAF_tumor'] / purity) * (purity * snv_cn['total_raw_e'] + 2 * (1 - purity))

    cnv['chromosome'] = cnv['chromosome'].replace({'x':'X', 'y':'Y'})

    fig = plt.figure(figsize=(15, 3))
    box = matplotlib.transforms.Bbox([[0.05, 0.05], [0.86, 0.95]])
    transform = matplotlib.transforms.BboxTransformTo(box)
    ax1 = fig.add_axes(transform.transform_bbox(box))

    major_col = 'major_raw'
    minor_col = 'minor_raw'
    ax1.set_ylabel('Raw copy number')
    chromosomes = [str(a) for a in list(range(1, 23)) + ['X', 'Y']]
    cn_plot.plot_cnv_genome(ax1, cnv, mincopies=-1, maxcopies=6, major_col=major_col, minor_col=minor_col, chromosomes=chromosomes, scatter=True)
    ax1.set_ylim((-0.5, None))
    ylim1 = ax1.get_ylim()
    # plot_scatter(snv_cn.coord, snv_cn.ccf, ax1, ylim=ylim1)
    plot_scatter(snv_cn.coord, snv_cn.snv_cn, ax1, ylim=ylim1) #
    label_colors = {'major_raw': 'tab:red', 'minor_raw':'tab:blue', 'SNV':'grey'}
    add_legend(fig, label_colors, marker='o')

    box = matplotlib.transforms.Bbox([[0.6, 0.05], [0.95, 0.95]])
    transform = matplotlib.transforms.BboxTransformTo(box)
    ax2 = fig.add_axes(transform.transform_bbox(box))
    # plot_hist(snv_cn.ccf, ax2, ylim=ylim1)
    plot_hist(snv_cn.snv_cn, ax2, ylim=ylim1) #
    
    fig.suptitle(sample)

    plt.tight_layout()
    break
    
```

```python
for sample in ['HCM-BROD-0199-C71-85A', 'HCM-BROD-0344-C25-85A', 'HCM-BROD-0195-C71-85A']:
    paths = meta[meta['isabl_sample_id']==sample]
    paths = paths[paths['result_type']=='command_log']
    assert paths.shape[0] == 1, paths
    path = paths['result_filepath'].iloc[0]
    out_dir, _ = os.path.split(path)
    inputs_path = f'{out_dir}/inputs.yaml'
    assert os.path.exists(inputs_path), inputs_path
    inputs = yaml.load(open(inputs_path))
    
    remixt_path = inputs['remixt']
    remixt = read_remixt.read(remixt_path)
    remixt['chrom'] = remixt['chrom'].str.replace('chr', '')
    remixt = remixt.dropna()
    
    somatic_calls_path = inputs['somatic_calls']
    somatic_calls = read_variant_calls.read(somatic_calls_path)
    somatic_calls['chrom'] = somatic_calls['chrom'].str.replace('chr', '')
    
    nsnvs = get_nygc_snvs(sample)
    somatic_calls = somatic_calls[somatic_calls.pos.isin(nsnvs.pos)]
    
    snv_copynumber = parse_snv_cn.parse(somatic_calls, remixt)
    
    cnv = remixt.copy()
    cnv.rename(columns={'chrom':'chromosome'}, inplace=True)
    cnv['chromosome_start'] = cnv.chromosome.str.upper().map(refgenome.info.chromosome_start)
    cnv['coord'] = cnv['start'] + cnv['chromosome_start']
    snv_cn = snv_copynumber.copy().dropna()
    snv_cn['chromosome_start'] = snv_cn.chrom.str.upper().map(refgenome.info.chromosome_start)
    snv_cn['coord'] = snv_cn['pos'] + snv_cn['chromosome_start']
    snv_cn['total_raw'] = snv_cn['major_raw'] + snv_cn['minor_raw']
    snv_cn['snv_cn'] = (snv_cn['VAF_tumor'] / purity) * (purity * snv_cn['total_raw'] + 2 * (1 - purity))
    snv_cn['snv_cn'] = snv_cn['snv_cn'].clip(upper=10)

    cnv['chromosome'] = cnv['chromosome'].replace({'x':'X', 'y':'Y'})

    fig = plt.figure(figsize=(15, 3))
    box = matplotlib.transforms.Bbox([[0.05, 0.05], [0.86, 0.95]])
    transform = matplotlib.transforms.BboxTransformTo(box)
    ax1 = fig.add_axes(transform.transform_bbox(box))

    major_col = 'major_raw'
    minor_col = 'minor_raw'
    ax1.set_ylabel('Raw copy number')
    chromosomes = [str(a) for a in list(range(1, 23)) + ['X', 'Y']]
    cn_plot.plot_cnv_genome(ax1, cnv, mincopies=-1, maxcopies=6, major_col=major_col, minor_col=minor_col, chromosomes=chromosomes, scatter=True)
    ax1.set_ylim((-0.5, None))
    ylim1 = ax1.get_ylim()
    # plot_scatter(snv_cn.coord, snv_cn.ccf, ax1, ylim=ylim1)
    plot_scatter(snv_cn.coord, snv_cn.snv_cn, ax1, ylim=ylim1) #
    label_colors = {'major_raw': 'tab:red', 'minor_raw':'tab:blue', 'SNV':'grey'}
    add_legend(fig, label_colors, marker='o')

    box = matplotlib.transforms.Bbox([[0.6, 0.05], [0.95, 0.95]])
    transform = matplotlib.transforms.BboxTransformTo(box)
    ax2 = fig.add_axes(transform.transform_bbox(box))
    # plot_hist(snv_cn.ccf, ax2, ylim=ylim1)
    plot_hist(snv_cn.snv_cn, ax2, ylim=ylim1) #
    
    fig.suptitle(sample)

    plt.tight_layout()
    # break
    
```

```python
snv_cn['snv_cn'].hist()
```

```python
fig, ax = plt.subplots(figsize=(4, 4))
tmp = snv_cn[snv_cn['chrom']=='7']
ax.scatter(x=tmp.coord, y=tmp.frac_cn)
ax.set_ylim((-1, 6));
```

```python
tmp
```

## QC of QC

```python
data.head()
```

```python

```

```python
data = snv_copynumber.copy()
```

```python
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
fig.suptitle(f'MSK: {sample}')
data['t_depth'].hist(bins=30, ax=ax1)
ax1.set_xlabel('t_depth')
ax1.set_xlim((0, 150))

data[
    (data['t_depth'] >= 10) 
    # &
    # (data['major_raw'].round(0) == 1) &
    # (data['minor_raw'].round(0) == 1)
]['VAF_tumor'].hist(bins=30, ax=ax2)
ax2.set_xlabel('depth>=10, all VAF')
ax2.set_xlim((0, 1))

data[
    (data['t_depth'] >= 20) 
    # &
    # (data['major_raw'].round(0) == 1) &
    # (data['minor_raw'].round(0) == 1)
]['VAF_tumor'].hist(bins=30, ax=ax3)
ax3.set_xlabel('depth>=20, all VAF')
ax3.set_xlim((0, 1))
```

### Comparison with NYGC

```python
vcf_path = '/juno/work/shah/users/chois7/tickets/hcmi/downloads/snvs/HCM-BROD-0199-C71-85A-01D-A786-36--HCM-BROD-0199-C71-10A-01D-A786-36.snv.indel.final.v6.annotated.vcf.gz'
vdf = pd.read_table(vcf_path, comment='#', names=['chrom', 'pos', 'ID', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'normal', 'tumor'])
vdf['chrom'] = vdf['chrom'].str.replace('chr', '')
```

```python
tdf = vdf['tumor'].str.split(':', expand=True)
tdf.columns = ['AD', 'DP', 'VAF']
tdf[['DP', 'VAF']] = tdf[['DP', 'VAF']].astype(float)
```

```python
sample
```

```python
data = tdf.copy()
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
fig.suptitle(f'NYGC: {sample}')
data['DP'].hist(bins=30, ax=ax1)
ax1.set_xlabel('t_depth')
ax1.set_xlim((0, 150))

data[
    (data['DP'] >= 10)
]['VAF'].hist(bins=30, ax=ax2)
ax2.set_xlabel('depth>=10, all VAF')
ax2.set_xlim((0, 1))

data[
    (data['DP'] >= 20)
]['VAF'].hist(bins=30, ax=ax3)
ax3.set_xlabel('depth>=20, all VAF')
ax3.set_xlim((0, 1))
```

```python
vdf
```

### position-matched

```python
msk = snv_copynumber.copy()
ngc = vdf.copy()
```

```python
both_ngc = ngc[
    ngc['pos'].isin(msk['pos'])
]

both_msk = msk[
    msk['pos'].isin(ngc['pos'])
]
```

```python
both_ngc.shape, both_msk.shape
```

```python
both_ngc
```

```python
data = both_msk.copy()

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
fig.suptitle(f'MSK: {sample}')
data['t_depth'].hist(bins=30, ax=ax1)
ax1.set_xlabel('t_depth')
ax1.set_xlim((0, 150))

data[
    (data['t_depth'] >= 10) 
    # &
    # (data['major_raw'].round(0) == 1) &
    # (data['minor_raw'].round(0) == 1)
]['VAF_tumor'].hist(bins=30, ax=ax2)
ax2.set_xlabel('depth>=10, all VAF')
ax2.set_xlim((0, 1))

data[
    (data['t_depth'] >= 20) 
    # &
    # (data['major_raw'].round(0) == 1) &
    # (data['minor_raw'].round(0) == 1)
]['VAF_tumor'].hist(bins=30, ax=ax3)
ax3.set_xlabel('depth>=20, all VAF')
ax3.set_xlim((0, 1))
```

```python
data = both_ngc.copy()
data = data['tumor'].str.split(':', expand=True)
data.columns = ['AD', 'DP', 'VAF']
data[['DP', 'VAF']] = data[['DP', 'VAF']].astype(float)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
fig.suptitle(f'NYGC: {sample}')
data['DP'].hist(bins=30, ax=ax1)
ax1.set_xlabel('t_depth')
ax1.set_xlim((0, 150))

data[
    (data['DP'] >= 10)
]['VAF'].hist(bins=30, ax=ax2)
ax2.set_xlabel('depth>=10, all VAF')
ax2.set_xlim((0, 1))

data[
    (data['DP'] >= 20)
]['VAF'].hist(bins=30, ax=ax3)
ax3.set_xlabel('depth>=20, all VAF')
ax3.set_xlim((0, 1))
```
