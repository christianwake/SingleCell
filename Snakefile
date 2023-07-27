
def maybe_skip_integration(qcdat):
  QC_files = defaultdict(str, zip(qcdat.step_name, qcdat.skip))
  if(QC_files['Integration'] == True):
    out = results_dir + 'Normalized.RDS'
  else:
    out = results_dir + 'PostQC4.RDS'
  return(out)

rule make_gtf_R_object:
  input:
    config['refgtf']
  params: 
    scripts = config['SC_scripts']
  output:
    'data/gtf.RDS'
  shell:
    'Rscript {params.scripts}/gtf.R {input} {output}'

rule QC_eval:
  input:
    "data/All_data.RDS"
    #QC_file ### So if the user edits this file to add or subtract an entire QC step, the DAG will be re-evaluating updating the following QC rules
  params:
    scripts = config['SC_scripts'],
    batch_candidates = config['batch'],
    test = config['test'],
    strat = config['stratifications'],
    plots_path = "plots/",
    results_dir = results_dir
  output:
    txt = results_dir + "batch_evaluation.txt",
    pdf = results_dir + "batch_evaluation.pdf",
    cp = touch('data/checkpoint.done') ### Makes rule 'Filter_whole_samples' dependent on this checkpoint
  shell:
    """
    mkdir -p {params.plots_path}
    mkdir -p {params.results_dir}
    Rscript {params.scripts}/QC_eval.R {input} '{params.batch_candidates}' '{params.test}' '{params.strat}' '{params.plots_path}' {output}
    touch {output.cp}
    """

### For the case that the later QC steps reveal some wholly bad samples
rule Filter_whole_samples:
  input:
    **maybe_skip_QC_stepN(qcdat, 1) ### Returns RDS file as 'rds' and file with QC details as 'qc'
  params:
    scripts = config['SC_scripts']
  output:
    results_dir + "PostQC1.RDS"
  shell:
    "Rscript {params.scripts}/Filter_samples.R {input.rds} {input.qc} {output}"

### MOD - optionally input gtf to better determine MT genes for quantification
rule Impute:
  input:
    **maybe_skip_QC_stepN(qcdat, 2) ### Returns RDS file as 'rds' and file with QC details as 'qc'
  params:
    scripts = config['SC_scripts']
  output:
    results_dir + "PostQC2.RDS",
    results_dir + "Imputed.pdf"
  resources: mem_mb=240000 
  shell:
    "Rscript {params.scripts}/Impute.R {input.rds} {input.qc} {output}"

rule Cell_and_feature_filter:
  input:
    **maybe_skip_QC_stepN(qcdat, 3), ### Returns RDS file as 'rds' and file with QC details as 'qc'
    gtf = 'data/gtf.RDS'
  output:
    results_dir + "PostQC3.RDS",
    results_dir + "Cell_filters.pdf",
    results_dir + "Filters_Nfeatures.txt",
    results_dir + "Excluded_genes.txt",
    results_dir + "Filters_Ncells.txt"
  params:
    scripts = config['SC_scripts']
  shell:
    "Rscript {params.scripts}/Cell_and_feature_filters.R {input.rds} {input.qc} '{input.gtf}' {output}"

rule Normalize:
  input: 
    **maybe_skip_QC_stepN(qcdat, 4, spec = False), ### Returns RDS file as 'rds' and file with QC details as 'qc'
    exclude = results_dir + "Excluded_genes.txt"
  params:
    scripts = config['SC_scripts'],
    norm = config['norm_method']
  output:
    results_dir + "Normalized.RDS",
  shell:
    "Rscript {params.scripts}/Normalization.R {input.rds} {input.exclude} '{params.norm}' {output}"

### Scale and Integrate:
rule Integrate:
  input:
    #**maybe_skip_QC_stepN(qcdat, 4), ### Returns RDS file as 'rds' and file with QC details as 'qc'
    rds = results_dir + "Normalized.RDS",
    qc = qcdat.iloc[4,0],
    exclude = results_dir + "Excluded_genes.txt"
  params:
    scripts = config['SC_scripts'],
    norm = config['norm_method']
  output:
    results_dir + "PostQC4.RDS",
    results_dir + "Integration.pdf"
  shell:
    "Rscript {params.scripts}/Integrate.R {input.rds} {input.qc} {input.exclude} {params.norm} {output}"

### Connects to the QC checkpoint
rule cluster:
  input:
    maybe_skip_integration(qcdat),
    exclude = results_dir + "Excluded_genes.txt",
    gtf = 'data/gtf.RDS'
  params:
    scripts = config['SC_scripts'],
    gene_file = config['gene_file'],
    clust_resolution = 'RNA-' + config['clust_resolution'],
    test = config['test'],
    integration = QC_specs['step4']
  output:
    results_dir + "Clustered.RDS",
    results_dir + "Clustered.pdf"
  resources: mem_mb=240000
  shell:
    "Rscript {params.scripts}/Clustering_SC.R {input} '{params.gene_file}' '{params.test}' '{params.integration}' {output} '{params.clust_resolution}'"

### Triggers DAG evaluation after successful execution of this rule
checkpoint cluster_filter:
  input:
    results_dir + "Clustered.RDS",
    QC_specs['step5'],
    'data/gtf.RDS'
  params:
    scripts = config['SC_scripts'],
  output:
    results_dir + "Cluster_filtered.RDS",
    results_dir + "Cluster_filters.txt",
    results_dir + "Cluster_filters.pdf"
  shell:
    "Rscript {params.scripts}/Cluster_filters.R {input} {output}"

### This rule may be skipped, after evaluation of cluster_filter checkpoint
rule post_filter_cluster:
  input:
    results_dir + "Cluster_filtered.RDS",
    results_dir + "Excluded_genes.txt",
    'data/gtf.RDS'
  params:
    scripts = config['SC_scripts'],
    gene_file = config['gene_file'],
    clust_resolution = 'RNA-' + config['clust_resolution2'],
    test = config['test'],
    integration = QC_specs['step4']
  output:
    results_dir + "Clustered_again.RDS",
    results_dir + "Clustered_again.pdf"
  shell:
    #"Rscript {params.scripts}/Clustering_SC.R {input} '{params.negative_markers}' '{params.gene_file}' '{params.clust_resolution}' '{params.test_variable}' '{params.integration}' {output}"
    "Rscript {params.scripts}/Clustering_SC.R {input} '{params.gene_file}' '{params.test}' '{params.integration}' {output} '{params.clust_resolution}'"

rule Cluster_FindAllMarkers:
  input:
    maybe_skip_second_cluster,
    results_dir + "Excluded_genes.txt",
    'data/gtf.RDS'
  params:
    scripts = config['SC_scripts'],
    results_dir = results_dir
  output:
    results_dir + "Cluster_DE/Cluster_DE.RDS",
    results_dir + "Cluster_DE/Cluster_DE.tsv"
  resources: mem_mb=240000
  shell:
    "Rscript {params.scripts}/Cluster_DE_SC.R {input} {output}"


rule fgsea_FindAllMarkers:
  input:
    maybe_skip_annotation(do_annot),
    results_dir + 'Cluster_DE/Cluster_DE.tsv',
    results_dir + "Excluded_genes.txt",
    'data/gtf.RDS'
  params: 
    scripts = config['SC_scripts'],
    gmt = config['pathways'],
    custom_sets = config['custom_sets'],
    species = config['species']
  output:
    results_dir + "Cluster_DE/Cluster_fgsea.tsv"
  shell:
    "Rscript {params.scripts}/fgsea.R {input} '{params.gmt}' '{params.species}' {output} '{params.custom_sets}'"

rule annotation:
  input:
    maybe_skip_second_cluster,
    results_dir + "Excluded_genes.txt",
    'data/gtf.RDS'
  params:
    scripts = config['SC_scripts'],
    annotation = config['annotation_reference'],
    method = config['annotation_method'],
    species = config['species'],
    cell_type = config['tissue']
  output:
    results_dir + "Mapped.RDS"
    #results_dir + "Mapping.pdf"
  resources: mem_mb=240000
  shell:
    "Rscript {params.scripts}/Annotation.R {input} '{params.annotation}' '{params.method}' '{params.species}' '{params.cell_type}' {output}"

rule FindAllMarkers_to_excel:
  input:
    results_dir + "Cluster_DE.RDS"
  params: 
    scripts = config['SC_scripts']
  output:
    results_dir + "Cluster_DE.xls"
  shell:
    "Rscript {params.scripts}/Excel_FindAllMarkers.R {input} {output} 0.05"

rule pseudobulk_DE:
  input:
    #results_dir + 'Cluster_filtered.RDS',
    maybe_skip_second_cluster,
    results_dir + "Excluded_genes.txt",
    'data/gtf.RDS'
  params: 
    scripts = config['SC_scripts'],
    path = results_dir + 'DE/{test_name}/'
  output:
    #results_dir + 'DE/{test}/{strat}/DE_results.tsv'
    results_dir + 'DE/{test_name}/{strat_name}/{test_name}-{value1}-{value2}_DE_within_{strat_name}-{strat_values}.tsv',
    results_dir + 'DE/{test_name}/{strat_name}/{test_name}-{value1}-{value2}_DE_within_{strat_name}-{strat_values}.pdf'
  shell:
    """
    mkdir -p {params.path}
    #Rscript {params.scripts}/pseudobulk_DE.R {input} {wildcards.test_name} {wildcards.strat_name} {output}
    Rscript {params.scripts}/pseudobulk_DE.R {input} '{wildcards.test_name}' '{wildcards.value1}' '{wildcards.value2}' '{wildcards.strat_name}' '{wildcards.strat_values}' {output}
    """

rule fgsea_test_DE:
  input:
    maybe_skip_second_cluster,
    #results_dir + 'DE/{test}/{strat}/DE_results.tsv',
    results_dir + 'DE/{test_name}/{strat_name}/{test_name}-{value1}-{value2}_DE_within_{strat_name}-{strat_values}.tsv',
    results_dir + "Excluded_genes.txt",
    'data/gtf.RDS'
  resources: mem_mb=240000
  params: 
    scripts = config['SC_scripts'],
    gmt = config['pathways'],
    custom_sets = config['custom_sets'],
    species = config['species']
  output:
    #results_dir + 'DE/{test}/{strat}/fgsea.tsv'
    results_dir + 'DE/{test_name}/{strat_name}/{test_name}-{value1}-{value2}_GSEA_within_{strat_name}-{strat_values}.tsv'
  shell:
    "Rscript {params.scripts}/fgsea.R {input} '{params.gmt}' '{params.species}' {output} '{params.custom_sets}'"

rule Excel_DE:
  input:
    [results_dir + "DE/{test_name}/{strat_name}/{test_name}-{value1}-{value2}_DE_within_{strat_name}-{strat_values}.tsv".format(test_name = comp_dict[dkey]['test_name'], value1 = comp_dict[dkey]['value1'], value2 = comp_dict[dkey]['value2'], strat_name = comp_dict[dkey]['strat_name'], strat_values = comp_dict[dkey]['strat_values']) for dkey in comp_dict.keys()]
  params:
    scripts = config['SC_scripts'],
    test = config['test']
  output:
    results_dir + "DE/DE.xls",
    results_dir + "DE/DE_sig.xls"
  shell:
    "Rscript {params.scripts}/Combine_results.R {output} 0.05 {input}"

#test1 = {'test_name': 'Celltype', 'val1':'CD4', 'val2':'CD8', 'strat_name':'All', 'strat_values':'All'}
#test2 = {'test_name': 'Time', 'val1':'30', 'val2':'31', 'strat_name':'Celltype', 'strat_values':'CD4'}
#comp_dict = {'test1': test1, 'test2': test2}
### Testing
rule Excel_fgsea:
  input:
    [results_dir + "DE/{test_name}/{strat_name}/{test_name}-{value1}-{value2}_GSEA_within_{strat_name}-{strat_values}.tsv".format(test_name = comp_dict[dkey]['test_name'], value1 = comp_dict[dkey]['value1'], value2 = comp_dict[dkey]['value2'], strat_name = comp_dict[dkey]['strat_name'], strat_values = comp_dict[dkey]['strat_values']) for dkey in comp_dict.keys()]
  params:
    scripts = config['SC_scripts'],
    test = config['test']
  output:
    results_dir + "DE/GSEA.xls",
    results_dir + "DE/GSEA_sig.xls"
  shell:
    "Rscript {params.scripts}/Combine_results.R {output} 0.05 {input}"

rule SubCluster:
  input:
    maybe_skip_annotation(do_annot),
    results_dir + "Excluded_genes.txt"
  params:
    scripts = config['SC_scripts'],
    results_dir = results_dir
  output:
    results_dir + "SubClusters/{cluster}.RDS"
  resources: mem_mb=240000
  shell:
    """
    mkdir -p {params.results_dir}/SubClusters/
    Rscript {params.scripts}/SubCluster.R {input} {wildcards.cluster} {output}
    """

rule Combine_SubClusters:
  input:
    expand(os.path.join(results_dir, 'SubClusters', "{cluster}.RDS"), cluster = sub_clusters)
  output:
    os.path.join(results_dir, "SubClusters", "SubClusters.csv")
  shell:
    "ls {input} > {output}"

rule App_data:
  input:
    ### maybe_skip_annotation(config['species'], config['tissue']),
    maybe_skip_annotation(config['annotation_reference']),
    results_dir + 'Cluster_DE.RDS',
    results_dir + "Excluded_genes.txt",
    'data/gtf.RDS',
    results_dir + "Cluster_DE/Singles.txt",
    results_dir + "Cluster_DE.RDS", ### FindAllMarkers
    #expand(os.path.join(results_dir, 'DE', '{test}', "{strat}", "DE_results.tsv"), test = tests, strat = strat_list), ### pseudobulk DE
    expand(os.path.join(results_dir, 'SubClusters', "{cluster}.RDS"), cluster = sub_clusters) ### 
  params: 
    scripts = config['SC_scripts'],
    project = config['project'],
    username = config['username'],
    server = config['server'],
    results_dir = results_dir,
    additional = config['additional_columns'],
    gene_file = config['gene_file']
  output:
    results_dir + "App/Data.RData",
    results_dir + "Extra_plots.pdf"
  shell:
    """
    mkdir -p {params.results_dir}/App/
    cp {params.scripts}/App.R {params.results_dir}/App/
    cp {params.scripts}/../../sc_functions.R {params.results_dir}/App/
    #Rscript {params.scripts}/Create_app_data.R {input} '{params.additional}' '{params.gene_file}' {output} '{params.project}' '{params.username}' '{params.server}'
    Rscript {params.scripts}/Create_app_data.R {output} '{params.additional}' '{params.gene_file}' '{params.project}' '{params.username}' '{params.server}' {input}
    """

