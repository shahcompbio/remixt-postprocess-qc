import os
import subprocess
import pandas as pd

if not os.path.exists(config['log_dir']): subprocess.run(f'mkdir -p {config["log_dir"]}', shell=True)
if not os.path.exists(config['tmp_dir']): subprocess.run(f'mkdir -p {config["tmp_dir"]}', shell=True)

#SAMPLES = ['HCM-CSHL-0058-C34-86A']
SAMPLES = pd.read_table(config['metadata'])['isabl_sample_id'].unique()

rule all:
    input:
        expand('results/remixt-pp-qc/{sample}.pdf', sample=SAMPLES),

def _get_remixtpp_cn_path(wildcards):
    tb_path = config['metadata']
    tb = pd.read_table(tb_path)
    tb = tb[tb['result_type']=='remixt_cn']
    tb = tb[tb['isabl_sample_id'] == wildcards.sample]
    assert tb.shape[0] == 1, f'tb=\n{tb}'
    path = tb['result_filepath'].iloc[0]
    return path

rule plot_qc:
    input:
        cn = _get_remixtpp_cn_path,
    output:
        pdf = 'results/remixt-pp-qc/{sample}.pdf',
    shell:
        'python scripts/plot_remixt_pp.py -i {input.cn} -o {output.pdf}'

