"""
This analysis searches for interesting genes/proteins that were differentially expressed in presence/absence of AC
But it also searches for proteins that may be involved in attachment/interactions with nanomaterials
"""

import pandas as pd
from paper_utils import parse_fasta

OUT = 'formicicum_proteomics'
SAMPLES = [f'CS{i}' for i in range(1, 7)]
FASTA_NAMES = [name.split('|')[1] for name in parse_fasta(f'{OUT}/formicicum_proteome.fasta').keys()]
up_res = pd.read_csv(f'{OUT}/results/upimapi/uniprotinfo.tsv', sep='\t')
recog_res = pd.read_excel(f'{OUT}/reCOGnizer_results.xlsx', sheet_name='COG')
recog_res = recog_res.groupby('qseqid')[recog_res.columns.tolist()[1:]].first().reset_index().rename(
    columns={'qseqid': 'Entry'})
recog_res['Entry'] = recog_res['Entry'].apply(lambda x: x.split('|')[1])

quant_df = pd.DataFrame(FASTA_NAMES, columns=['Entry'])
for sample in SAMPLES:
    ps_report = pd.read_csv(
        f'{OUT}/{sample}/experiment_{sample}_1_Default_Protein_Report.txt', sep='\t', index_col=0,
        usecols=['Main Accession', '#Validated PSMs']).reset_index().rename(
        columns={'#Validated PSMs': sample, 'Main Accession': 'Entry'})
    quant_df = pd.merge(quant_df, ps_report, how='left', on='Entry')

all_info = pd.merge(up_res, recog_res, on='Entry', how='left')
all_info = pd.merge(all_info, quant_df, on='Entry', how='left')
all_info[SAMPLES] = all_info[SAMPLES].fillna(value=0.0)

