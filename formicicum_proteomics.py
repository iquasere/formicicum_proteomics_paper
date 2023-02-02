import pandas as pd
import os
import numpy as np
import pathlib


database = 'task1/database.fasta'
decoy_database = database.replace('.fasta', '_concatenated_target_decoy.fasta')
if not os.path.isfile(decoy_database):
    run_command('searchgui eu.isas.searchgui.cmd.FastaCLI -in {} -decoy'.format(database))
else:
    print(decoy_database + ' already exists!')

def generate_parameters_file(output, database):
    run_pipe_command(('searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI -out {} -db {} -prec_tol 10 '
                '-frag_tol 0.02 -enzyme Trypsin -fixed_mods "Carbamidomethylation of C" -variable_mods "Oxidation '
                'of M, Acetylation of protein N-term" -mc 2').format(output, database))
                
def peptide_spectrum_matching(spectra_folder, output, parameters_file, threads='12',
                              search_engines = ['xtandem', 'myrimatch', 'msgf'], max_memory=4096):
    run_command('searchgui eu.isas.searchgui.cmd.SearchCLI -Xmx{}M -spectrum_files {} -output_folder {} -id_params {} -threads {}{}'.format(
        max_memory, spectra_folder, output, parameters_file, threads,
        ''.join([' -{} 1'.format(engine) for engine in search_engines])))
        
def browse_identification_results(spectra_folder, parameters_file,
                                  searchcli_output, peptideshaker_output, experiment_name='experiment',
                                  sample_name='sample', replicate_number='1', max_memory=4096):
    try:
        run_command(('peptide-shaker -Xmx{} eu.isas.peptideshaker.cmd.PeptideShakerCLI -spectrum_files {} ' +
                   '-experiment {} -sample {} -replicate {} -identification_files {} -out {}').format(
        max_memory, spectra_folder, experiment_name, sample_name, replicate_number, searchcli_output,
        peptideshaker_output))
    except:
        print('Producing Peptide-Shaker result failed! Maybe no identifications were obtained?')
        
def generate_reports(peptideshaker_output, reports_folder, reports_list=[str(n) for n in range(12)]):
    print('Created ' + reports_folder)
    pathlib.Path(reports_folder).mkdir(parents=True, exist_ok=True)                 # creates folder for reports
    run_command('peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in {} -out_reports {} -reports {}'.format(
        peptideshaker_output, reports_folder, ','.join(reports_list)))

def run_compomics(number):
    folder = 'task1/CS{}'.format(letter)
    generate_parameters_file('params.par', decoy_database)
    peptide_spectrum_matching(folder, folder, 'params.par')
    browse_identification_results(folder, 'params.par', folder + '/searchgui_out.zip', folder + '/ps_output.cpsx', sample_name='CS' + number)
    
for i in range(1,7):
    run_compomics(i)

def get_cs(number):
    res = pd.read_csv('task1/CS{0}/experiment_CS{0}_1_Default_Protein_Report_with_non-validated_matches.txt'.format(number),sep='\t')
    res = res[['Main Accession', '#Validated PSMs']]
    res.columns = ['Main Accession','CS{} (#Validated PSMs)'.format(number)]
    return res

cs1 = get_cs(1)
cs2 = get_cs(2)
cs3 = get_cs(3)
cs4 = get_cs(4)
cs5 = get_cs(5)
cs6 = get_cs(6)

css = pd.DataFrame(columns = ['Main Accession'])
for cs in [cs1,cs2,cs3,cs4,cs5,cs6]:
    css = pd.merge(css, cs, on='Main Accession',how='outer')
css.fillna(value=0.0,inplace=True)
css[['CS{} (#Validated PSMs)'.format(number) for number in [1,2,3,4,5,6]]] = css[['CS{} (#Validated PSMs)'.format(number) for number in [1,2,3,4,5,6]]].astype(int)


all = pd.DataFrame(columns=['qseqid','sseqid','evalue','DB ID','DB description','Sequence','product_name','ec_number'])
for sheet_name in ['Pfam','NCBIfam','Protein_Clusters','Smart','TIGRFAM']:#'COG','KOG','CDD',
    print(sheet_name)
    data = pd.read_excel('reCOGnizer_results_top15.xlsx',sheet_name=sheet_name)
    #data = data[data['evalue']<=0.01]
    if sheet_name == 'COG':
        data = data[data['COG protein description'].notnull()]
        data = data[data['DB description'].str.contains('quorum') | data['COG protein description'].str.contains('quorum')][['qseqid','sseqid','evalue','DB ID','DB description','Sequence','COG protein description','EC number']]
        data.columns = data.columns.tolist()[:-2] + ['product_name','ec_number']
    elif sheet_name == 'KOG':
        data = data[data['KOG protein description'].notnull()]
        data = data[(data['DB description'].str.contains('quorum')) | (data['KOG protein description'].str.contains('quorum'))][['qseqid','sseqid','evalue','DB ID','DB description','Sequence','KOG protein description']]
        data.columns = data.columns.tolist()[:-1] + ['product_name']
        data['ec_number'] = [np.nan] * len(data)
    elif sheet_name == 'Smart':
        data = data[data['Smart description'].notnull()]
        data = data[(data['DB description'].str.contains('quorum')) | (data['Smart description'].str.contains('quorum'))][['qseqid','sseqid','evalue','DB ID','DB description','Sequence','Smart description']]
        data.columns = data.columns.tolist()[:-1] + ['product_name']
        data['ec_number'] = [np.nan] * len(data)
    else:
        data['product_name'].fillna(value='',inplace=True)
        data = data[data['DB description'].notnull()]
        data = data[(data['DB description'].str.contains('quorum')) | (data['product_name'].str.contains('quorum'))][['qseqid','sseqid','evalue','DB ID','DB description','Sequence','product_name','ec_number']]
    all = pd.concat([all, data])
all.to_excel('qs.xlsx')
    
qs = pd.read_excel('qs.xlsx', sheet_name='Confident QS')
qs['Main Accession'] = [ide.split('|')[1] for ide in qs.qseqid]
result = pd.merge(qs, css, on='Main Accession', how='left')
for i in range(1,7):
    qs['CS{} (normalized)'.format(i)]=qs['CS{} (#Validated PSMs)'.format(i)] / qs['CS{} (#Validated PSMs)'.format(i)].sum() * 10000

# wget https://github.com/qhmu/QSDB/raw/master/data/QSDB.aa.fasta.gz -P resources_directory/
# gunzip resources_directory/QSDB.aa.fasta.gz


