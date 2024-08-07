import pandas as pd
import subprocess
import os
import re
from collections import defaultdict
import csv
import sys

#project_config_file = sys.argv[1]
#print(project_config_file)
configfile: 'project_config_file.yaml'
#configfile: project_config_file
### Read configuration file
config = defaultdict(str, config)

print(config['QC_name'])
QC_name = config['QC_name']
results_dir = os.path.join('results', QC_name) + '/'
log_dir = os.path.join(results_dir, 'logs/')

### Should we attempt automated cell type annotation?
do_annot = False
if(config['annotation_reference'] != ''):
  do_annot = True

### Which data types are included?
assays_list = ['RNA']
clusters_list = ['RNA_clusters']
if(config['Data_type'] == 'CITESeq'):
  assays_list = assays_list + ['prot']
  clusters_list = clusters_list + ['prot_clusters', 'wnn_clusters']
if(do_annot):
  clusters_list = clusters_list + ['predicted.celltype.l1']

### Is the data sample-hashed?
dehash = config['hashed']
if dehash:
  methods = config['dehash_method_comparisons'].split(',')
  dehash_calls_tsv = defaultdict(str)
  dehash_CRint_calls_tsv = defaultdict(str)
  dehash_threshs_csv = defaultdict(str)
  for method in methods:
    ### If the input method is a file, lets assume they are hashing results from some custom method. For now, we have to assume the cell IDs will match those in the Cell Ranger outputs
    if(os.path.isfile(method)):
      dehash_threshs_csv[method] = method
      call_file = results_dir + 'Dehash/' + os.path.basename(re.sub("Dehash_threshs_", 'Dehash_calls_', re.sub('csv$', 'tsv', method)))
      dehash_calls_tsv[method] = call_file
      call_file = results_dir + 'Dehash/' + os.path.basename(re.sub("Dehash_threshs_", 'Dehash_CRint_calls_', re.sub('csv$', 'tsv', method)))
      dehash_CRint_calls_tsv[method] = call_file
    else: ### Otherwise, define the result output based on the input method name
      dehash_calls_tsv[method] = results_dir + 'Dehash/Dehash_calls_' + method + '.tsv'
      dehash_CRint_calls_tsv[method] = results_dir + 'Dehash/Dehash_CRint_calls_' + method + '.tsv'
      dehash_threshs_csv[method] = results_dir + 'Dehash/Dehash_threshs_' + method + '.csv'
else:
  dehash_calls_tsv = {}
  dehash_threshs_csv = {}
  dehash_CRint_calls_tsv = {}

print(dehash_CRint_calls_tsv)

### If not input, assume the method to continue with is the first one input in the string of options
if(config['dehash_method'] == ''):
  utilized_method_name = config['dehash_method_comparisons'].split(',')[0]
else:
  utilized_method_name = config['dehash_method']
print(utilized_method_name)


### If the second cluster resolution is not specified in the config file, use the first cluster resolution.
if 'RNA_resolution2' not in config.keys():
  config['RNA_resolution2'] = config['RNA_resolution']
if 'prot_resolution2' not in config.keys():
  config['prot_resolution2'] = config['prot_resolution']

comp_dict = defaultdict(int)
with open(config['comparison_file'], newline = '') as csvfile:
  reader = csv.DictReader(csvfile)
  keyint = 0
  for row in reader:
    comp_dict[keyint] = row
    ### If no input strat_name, set dictionary values to 'All'
    if(comp_dict[keyint]['strat_name1'] == ''):
      comp_dict[keyint]['strat_name1'] = 'All'
      comp_dict[keyint]['strat_values1A'] = 'All'
      comp_dict[keyint]['strat_values1B'] = 'All'
    else: ### If strat_name is input, and strat_valueA is input but strat_value B is not, set the second value to the first
      if(comp_dict[keyint]['strat_values1A'] != '' and comp_dict[keyint]['strat_values1B'] == ''):
        comp_dict[keyint]['strat_values1B'] = comp_dict[keyint]['strat_values1A']
    ### If no input strat_name, set dictionary values to 'All'
    if(comp_dict[keyint]['strat_name2'] == ''):
      comp_dict[keyint]['strat_name2'] = 'All'
      comp_dict[keyint]['strat_values2A'] = 'All'
      comp_dict[keyint]['strat_values2B'] = 'All'
    else: ### If strat_name is input, and strat_valueA is input but strat_value B is not, set the second value to the first
      if(comp_dict[keyint]['strat_values2A'] != '' and comp_dict[keyint]['strat_values2B'] == ''):
        comp_dict[keyint]['strat_values2B'] = comp_dict[keyint]['strat_values2A']
    ### If either test value is empty, set to NA
    ### If either test value is empty, set to NA
    if(comp_dict[keyint]['value1'] == ''): 
      comp_dict[keyint]['value1'] = 'NA'
    if(comp_dict[keyint]['value2'] == ''): 
      comp_dict[keyint]['value2'] = 'NA'
    keyint += 1

### Read list of clusters to sub-cluster
sc_file = os.path.join(os.getcwd(), config['sub_cluster'])
if os.path.exists(sc_file):
  sc = pd.read_csv(sc_file)
  ### Because cluster names are likely numbers only
  sc['cluster'] = sc['cluster'].astype(str)
sub_clusters = ['-'.join(list(row)) for index,row in sc.iterrows()]
print('Sub clustering: ')
print(sub_clusters)

### Read QC summary file (input to batch_eval checkpoint) and creates dictionary to hold the file paths held within
QC_file = os.path.join(os.getcwd(), config['QC_file'])
qcdat = pd.read_csv(QC_file)
#### Add missing steps to qcdat with '' file column, including step 0 (no QC done yet)
#step_names = ['Base data', 'DSB background', 'Sample filter', 'Imputation', 'Downsample', 'Transcript filter', 'Cell filter', 'Integration', 'Feature filter', 'Cluster filter']
#step_names = ['Base data', 'Sample filter', 'Imputation', 'Cell filter', 'Integration', 'Cluster filter']
steps = list(set(list(range(len(step_names)))) - set(qcdat.step))
d = {'step':steps, 'file':['' for s in steps]}
### If there are steps to add, add them. Conditional is necessary because merging an empty dictionary converts step numbers from str to num (for unknown reason)
if(len(steps) > 0):
  qcdat = pd.concat([qcdat, pd.DataFrame(d)])
### Reorder rows
qcdat['step_num'] = ['step' + str(q) for q in qcdat.step]
qcdat['post_name'] = [results_dir + 'PostQC' + str(q) + '.RDS' for q in qcdat.step]
qcdat = qcdat.set_index('step')
qcdat = qcdat.sort_index(ascending = True)
### Add step_name
qcdat['step_name'] = step_names
### Add original RDS file
#qcdat.iloc[0, qcdat.columns.get_loc('post_name')] = 'data/All_data.RDS'
qcdat.iloc[0, qcdat.columns.get_loc('post_name')] = results_dir + 'All_data.RDS'
### Add skip information
#qcdat['skip'] = [f == '' for f in qcdat.file]
qcdat['skip'] = [not os.path.exists(f) for f in qcdat.file]
qcdat.iloc[0, qcdat.columns.get_loc('skip')] = False

print(qcdat)
QC_specs = defaultdict(str, zip(qcdat.step_num, qcdat.file))

### For RStudio Connect publication ...sigh
if('username' not in config.keys()): config['username'] = 'NA'
if('server' not in config.keys()): config['server'] = 'NA'

### For dynamic resource requests
size_step = [8000, 15000, 64000, 240000]
def get_mem_mb(wildcards, attempt):
  return size_step[attempt]

### If the last line of the second output of cluster_filters (Filters.txt) specifies that 0 cells were filtered, skip re-clustering
def maybe_skip_second_cluster(wildcards):
  checkp = checkpoints.cluster_filter.get().output[1]
  print(checkp)
  with open(checkp) as f:
    for line in f:
      if(re.search('Total removed', line)):
        filtered_line = line
  print(filtered_line)
  n_filtered = filtered_line.split('\t')[-1] ### Last item in the line is n
  if(int(n_filtered) == 0):
    print('Skip second cluster.')
    return results_dir + 'Cluster_filtered.RDS'
  else:
    print('Do second cluster.')
    return results_dir + 'Clustered_again.RDS'

def maybe_skip_QC_stepN(qcdat, stepN, spec = True):
  QC_files = defaultdict(str, zip(qcdat.step_num, qcdat.file))
  ### Name of the output RDS file for this QC step
  QC_rds = dict(zip(qcdat.step_num, qcdat.post_name))
  ### File with specifics for this step
  specs_file = QC_files['step' + str(stepN)]
  ### Slice the data frame to include only up to input stepN
  qcdat = qcdat[:stepN]
  ### Select the step name of the last step that is not skipped (column 'skip' = F). Use that step name as key in the RDS dictionary.
  rds_file = QC_rds[qcdat[qcdat.skip == False].iloc[-1, qcdat.columns.get_loc('step_num')]]
  out = {'rds':rds_file, 'touch':'data/checkpoint.done'}
  if spec:
    out['qc'] = specs_file
  return(out)

def maybe_skip_annotation(do_annot):
  if(do_annot):
    out = results_dir + 'Mapped.RDS'
    print('Doing annotation')
  else:
    out = maybe_skip_second_cluster
  return out

def maybe_skip_integration(qcdat):
  QC_files = defaultdict(str, zip(qcdat.step_name, qcdat.skip))
  if(QC_files['Integration'] == True):
    out = results_dir + 'RNA_normalized.RDS'
  else:
    out = results_dir + 'PostQC4.RDS'
  return(out)


