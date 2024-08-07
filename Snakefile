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
    results_dir + "All_data.RDS"
    #QC_file ### So if the user edits this file to add or subtract an entire QC step, the DAG will be re-evaluating updating the following QC rules
  params:
    scripts = config['SC_scripts'],
    test_file = config['comparison_file'],
    batch_candidates = config['batch'],
    plots_path = "plots/",
    results_dir = results_dir
  output:
    txt = results_dir + "batch_evaluation.txt",
    pdf = results_dir + "batch_evaluation.pdf",
    cp = touch('data/checkpoint.done') ### Makes rule 'Filter_whole_samples' dependent on this checkpoint
  resources: mem_mb=get_mem_mb
  shell:
    """
    mkdir -p {params.plots_path}
    mkdir -p {params.results_dir}
    Rscript {params.scripts}/QC_eval.R {input} '{params.test_file}' '{params.batch_candidates}' '{params.plots_path}' {output}
    touch {output.cp}
    """

### QC step0 (not optional for CITESeq)
### Background cells are those that are called that way both by CellRanger and by the dehashing step, and then are filtered further.
rule Create_DSB_background:
  input:
    samples = sample_sheet,
    tsv = dehash_calls_tsv[utilized_method_name],
    qc = QC_specs['step0'], 
    gtf = 'data/gtf.RDS'
  params:
    scripts = config['CITESeq_scripts'],
    runs_dir = config['runs_dir'],
    project = config['project'],
    results_dir = results_dir
  output:
    results_dir + "Background.RDS",
    results_dir + "Background.pdf"
  shell:
    """
    mkdir -p {params.results_dir}
    Rscript {params.scripts}/Create_background.R '{params.runs_dir}' '{params.project}' {input} {output}
    """

### QC step1, sample filter e.g. entire batches gone wrong
rule Filter_whole_samples_for_Seurat:
  input:
    **maybe_skip_QC_stepN(qcdat, 1) ### Returns RDS file as 'rds' and file with QC details as 'qc', and 'touch'
    #"data/All_data.RDS"
  params:
    scripts = config['SC_scripts'],
    samples = samples
  output:
    rds = results_dir + "PostQC1.RDS"
    #csv = results_dir + "Sample_sheet.csv"
  resources: mem_mb=240000
  run:
    shell('Rscript {params.scripts}/Filter_samples.R {input.rds} {input.qc} {output.rds}')

rule Filter_whole_samples_for_Sample_sheet: 
  input:
    sample_sheet = sample_sheet 
  output:
    csv = results_dir + "Sample_sheet.csv"
  run:
    samples = pd.read_table(input['sample_sheet'], sep = ',').set_index('Sample_ID', drop = False)
    qc_skip = qcdat[qcdat.step_num == 'step1'].skip.values[0]
    ### If Sample filter file is input, subset and print sample sheet into the results/qc_name directory
    if(qc_skip == False):
      print('Filtering samples')
      qc_file = qcdat[qcdat.step_num == 'step1'].file.values[0]
      ### subset samples
      dat = pd.read_table(qc_file, sep = ',')
      dat = dat[dat['type'] == 'batch_remove']
      og_size = samples.size
      ### For each line in the input csv file
      for index, row in dat.iterrows():
        ### Determine whether the line defines items to include or exlucde
        exclude = '!' in row['value']
        val_str = re.sub('!', '', row['value'])
        ### If the lines values include a comma, it will be considered a list of values
        if(',' in val_str):
          val_list = val_str.split(',')
        else:
          val_list = [val_str]
        for v in val_list:
          mid_size = samples.size      
          if(exclude):
            samples = samples[samples[row['feature']] != v]
          else: 
            samples = samples[samples[row['feature']] == v]
          print(str(mid_size - samples.size) + ' matching ' + v)
      print('Total of ' + str(og_size - samples.size) + ' removed')
    ### Write to output file
    samples.to_csv(output['csv'])

checkpoint set_batches:
  input:
    sample_sheet = results_dir + "Sample_sheet.csv"
  output:
    batch_dir = directory(results_dir + 'batches/')
  run:
    ### Create directory
    os.mkdir(output['batch_dir'])
    samples = pd.read_table(input['sample_sheet'], sep = ',').set_index('Sample_ID', drop = False)
    ### The batches based on the sample sheet output from the previous step (filtered)
    batchs = list(set(samples[config['batch']]))
    for batch in batchs:
      ### Create sub directory
      os.mkdir(os.path.join(output['batch_dir'], batch))
      ### Get sample sheet subset
      sub = samples[samples[config['batch']] == batch]
      ### Print
      sub.to_csv(os.path.join(output['batch_dir'], batch, 'Sample_sheet.csv'))

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

### QC step 3 to downsample
rule Downsample:
  input:
    **maybe_skip_QC_stepN(qcdat, 3), ### Returns RDS file as 'rds' and file with QC details as 'qc', and 'touch'
  params:
    scripts = config['CITESeq_scripts'],
    results_dir = results_dir
  output:
    results_dir + "PostQC3.RDS",
  resources: mem_mb=get_mem_mb
  shell:
    """
    Rscript {params.scripts}/Downsample.R {input.rds} {input.qc} {output}
    """

### Add output with N genes per filter criteria
rule Feature_filters:
  input:
    **maybe_skip_QC_stepN(qcdat, 4),
  params:
    scripts = config['SC_scripts'],
    results_dir = results_dir
  output:
    results_dir + "Excluded_genes.txt",
    results_dir + "Excluded_genes_analysis_only.txt",
    results_dir + "PostQC4.RDS",
    results_dir + "Excluded_genes.tsv"
  resources: mem_mb=get_mem_mb
  log:
    log_dir + "Feature_filters.log"
  shell:
    """
    Rscript {params.scripts}/Transcript_filters.R {input.rds} {input.qc} {output} > {log}
    """

### QC step 5, remove cells based on RNA content
rule Cell_filters:
  input:
    **maybe_skip_QC_stepN(qcdat, 5), ### Returns RDS file as 'rds' and file with QC details as 'qc', and 'touch'
    samples = results_dir + "Sample_sheet.csv",
    exclude = results_dir + "Excluded_genes.txt",
    gtf = 'data/gtf.RDS'
  params:
    scripts = config['SC_scripts'],
    batch_name = config['batch'],
    results_dir = results_dir
  output:
    results_dir + "batches/{batch}/Cell_filtered.RDS",
    results_dir + "batches/{batch}/Cell_filtered.pdf",
    results_dir + "batches/{batch}/Cell_filtered.txt"
  resources: mem_mb=get_mem_mb
  log:
    results_dir + "batches/{batch}/Cell_filtered.log"
  shell:
    """
    mkdir -p {params.results_dir}/batches/{wildcards.batch}/
    Rscript {params.scripts}/Cell_filters.R {input.rds} {input.qc} {input.samples} {input.exclude} {input.gtf} '{wildcards.batch}' '{params.batch_name}' {output}
    """

### Per batch, run DSB protein normalization
rule DSB:
  input:
    #runs_dir = config['runs_dir'],
    sdat = results_dir + "batches/{batch}/Cell_filtered.RDS",
    background = results_dir + "Background.RDS",
    samples = results_dir + "Sample_sheet.csv",
    #samples = sample_sheet,
    #csv = "data/Cell_data.csv",
    tsv = dehash_calls_tsv[utilized_method_name],
    #"data/Unnormalized_data.RDS"
  output:
    results_dir + "batches/{batch}/DSB_normalized_data.RDS"
  params:
    scripts = config['CITESeq_scripts'],
    batch_name = config['batch'],
    #htos = config['hashtags'],
    isocont = config['isotype_controls']
  shell:
    "Rscript {params.scripts}/DSB.R {input} '{wildcards.batch}' '{params.batch_name}' '{params.isocont}' {output}"

rule Combine_batches:
  input:
    aggregate_rds_inputs
    #[os.path.join(results_dir, str(Batch), 'RNA_cell_filtered.RDS') for Batch in batchs],
    #[os.path.join(results_dir, str(Batch), 'DSB_normalized_data.RDS') for Batch in batchs]
  params:
    scripts = config['CITESeq_scripts'],
    batch_name = config['batch'],
    isocont = config['isotype_controls']
  output:
    results_dir + "Combined_batches.RDS",
    results_dir + "Combined_batches.pdf"
  resources: mem_mb=240000
  log:
    log_dir + "Combine_batches.log"
  shell:
    """
    Rscript {params.scripts}/DSB_combine.R {output} {params.batch_name} '{params.isocont}' {input}
    """

rule RNA_Normalization:
  input: 
    #**maybe_skip_QC_stepN(qcdat, 4, spec = False), ### Returns RDS file as 'rds' and file with QC details as 'qc'
    rds = ancient(results_dir + "Combined_batches.RDS"),
    exclude = results_dir + "Excluded_genes.txt"
  params:
    scripts = config['SC_scripts'],
    norm = config['norm_method']
  output:
    results_dir + "RNA_normalized.RDS",
  resources: mem_mb=240000
  shell:
    "Rscript {params.scripts}/RNA_normalization.R {input.rds} {input.exclude} '{params.norm}' {output}"

rule batch_evaluation_again:
  input:
    #results_dir + "Combined_batches.RDS"
    results_dir + "RNA_normalized.RDS"
  params:
    scripts = config['SC_scripts'],
    test_file = config['comparison_file'],
    batch = config['batch'],
    plots_path = results_dir + '/post_DSB/'
  output:
    txt = results_dir + 'Batch.txt',
    pdf1 = results_dir + 'Post_QC.pdf'
    #pdf2 = results_dir + 'Batch.pdf'
  resources: mem_mb=get_mem_mb
  shell:
    """
    mkdir -p {params.plots_path}
    Rscript {params.scripts}/QC_eval.R {input} '{params.test_file}' '{params.batch}' '{params.plots_path}' {output}
    """

rule Combine_QC_outs:
  input:
    aggregate_txt_inputs
    #[os.path.join(results_dir, str(Batch), 'Cell_filtered.txt') for Batch in batchs],
    #[os.path.join(results_dir,str(Batch), 'Excluded_genes.txt') for Batch in batchs]
  params:
    scripts = config['SC_scripts'],
  output:
    results_dir + 'Cell_filtered.txt',
    results_dir + 'Batches_cells_remaining.txt',
    results_dir + 'Cell_filtered.xls',
    results_dir + 'Cell_filtered.pdf'
  shell:
    "Rscript {params.scripts}/Combine_QC_outputs.R {output} {input}"

### Scale and Integrate:
rule Integrate:
  input:
    **maybe_skip_QC_stepN(qcdat, 6), ### Returns RDS file as 'rds' and file with QC details as 'qc'
    exclude = results_dir + "Excluded_genes.txt"
  params:
    scripts = config['SC_scripts'],
    norm = config['norm_method']
  output:
    results_dir + "PostQC6.RDS",
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
    RNA_resolution = 'RNA-' + str(config['RNA_resolution']),
    prot_resolution = 'prot-' + str(config['prot_resolution']),
    test = config['comparison_file'],
    integration = QC_specs['step6']
  output:
    results_dir + "Clustered.RDS",
    results_dir + "Clustered.pdf"
  resources: mem_mb=240000
  shell:
    "Rscript {params.scripts}/Clustering.R {input} '{params.gene_file}' '{params.test}' '{params.integration}' {output} '{params.RNA_resolution}' '{params.prot_resolution}'"

### Triggers DAG evaluation after successful execution of this rule
checkpoint cluster_filter:
  input:
    results_dir + "Clustered.RDS",
    QC_specs['step7'],
    'data/gtf.RDS'
  params:
    scripts = config['SC_scripts'],
  output:
    results_dir + "Cluster_filtered.RDS",
    results_dir + "Cluster_filtered.txt",
    results_dir + "Cluster_filtered.pdf"
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
    RNA_resolution = 'RNA-' + config['clust_resolution2'],
    prot_resolution = 'prot-' + str(config['prot_resolution2']),
    test = config['comparison_file'],
    integration = QC_specs['step6']
  output:
    results_dir + "Clustered_again.RDS",
    results_dir + "Clustered_again.pdf"
  shell:
    #"Rscript {params.scripts}/Clustering_SC.R {input} '{params.negative_markers}' '{params.gene_file}' '{params.clust_resolution}' '{params.test_variable}' '{params.integration}' {output}"
    "Rscript {params.scripts}/Clustering.R {input} '{params.gene_file}' '{params.test}' '{params.integration}' {output} '{params.RNA_resolution}' '{params.prot_resolution}'"

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
    #results_dir + "Mapped.pdf"
  resources: mem_mb=240000
  shell:
    "Rscript {params.scripts}/Annotation.R {input} '{params.annotation}' '{params.method}' '{params.species}' '{params.cell_type}' {output}"

### FindAllMarkers Seurat function for cluster DE (one v all)
rule FAM_DE:
  input:
    maybe_skip_annotation(do_annot),
    results_dir + "Excluded_genes.txt",
    results_dir + "Excluded_genes_analysis_only.txt",
    'data/gtf.RDS'
  params:
    scripts = config['SC_scripts'],
    results_dir = results_dir
  output:
    results_dir + "Cluster_DE/DE_{assay}-{cluster_type}.RDS",
    results_dir + "Cluster_DE/DE_{assay}-{cluster_type}.tsv"
  resources: mem_mb=240000
  shell:
    """
    mkdir -p {params.results_dir}/Cluster_DE/
    Rscript {params.scripts}/Cluster_DE_SC.R {input} {output} {wildcards.assay} {wildcards.cluster_type}
    """

#rule FAM_GSEA:
#  input:
#    maybe_skip_annotation(do_annot),
#    results_dir + 'Cluster_DE/DE_{assay}-{cluster}.tsv',
#    results_dir + "Excluded_genes.txt",
#    'data/gtf.RDS'
#  params: 
#    scripts = config['SC_scripts'],
#    gmt = config['pathways'],
#    custom_sets = config['custom_sets'],
#    species = config['species']
#  output:
#    results_dir + "Cluster_DE/GSEA_{assay}-{cluster}.tsv",
#    results_dir + "Cluster_DE/GSEA_{assay}-{cluster}.pdf"
#  shell:
#    "Rscript {params.scripts}/fgsea.R {input} '{params.gmt}' '{params.species}' {output} '{params.custom_sets}'"

rule FAM_DE_to_excel:
  input:
    results_dir + "Cluster_DE/DE_{assay}-{cluster_type}.RDS"
  params: 
    scripts = config['SC_scripts']
  output:
    results_dir + "Cluster_DE/DE_{assay}-{cluster_type}.xls"
  shell:
    "Rscript {params.scripts}/FAM_DE_to_excel.R {input} {output} 0.05"

#rule FAM_GSEA_to_excel: 
#  input:
#    results_dir + "Cluster_DE/GSEA_{assay}-{cluster}.RDS"
#  params: 
#    scripts = config['SC_scripts']
#  output:
#    results_dir + "Cluster_DE/GSEA_{assay}-{cluster}.xls"
#  shell:
#    "Rscript {params.scripts}/FAM_GSEA_to_excel.R {input} {output} 0.05"

rule FAM_combine:
  input:
    expand(os.path.join(results_dir, 'Cluster_DE', 'DE_{assay}-{cluster_type}.RDS'), assay = assays_list, cluster_type = clusters_list)
    #expand(os.path.join(results_dir, 'Cluster_DE', 'GSEA_{assay}-{cluster_type}.RDS'), assay = assays_list, cluster_type = clusters_list)
  output:
    results_dir + 'Cluster_DE/FAM_files.txt'
  shell:
    "ls {input} > {output}"

rule pseudobulk_DE:
  input:
    maybe_skip_annotation(do_annot),
    results_dir + "Excluded_genes.txt",
    results_dir + "Excluded_genes_analysis_only.txt",
    'data/gtf.RDS'
  params: 
    scripts = config['SC_scripts'],
    path = results_dir + 'DE/{test_name}/{strat_name1}/'
  output:
    results_dir + 'DE/{test_name}/{strat_name1}/{test_name}-{value1}-{value2}_DE_Strat1-{strat_name1}-{strat_values1A}-{strat_values1B}_Strat2-{strat_name2}-{strat_values2A}-{strat_values2B}.tsv',
    results_dir + 'DE/{test_name}/{strat_name1}/{test_name}-{value1}-{value2}_DE_Strat1-{strat_name1}-{strat_values1A}-{strat_values1B}_Strat2-{strat_name2}-{strat_values2A}-{strat_values2B}.pdf'
  shell:
    """
    mkdir -p {params.path}
    Rscript {params.scripts}/pseudobulk_DE.R {input} '{wildcards.test_name}' '{wildcards.value1}' '{wildcards.value2}' '{wildcards.strat_name1}' '{wildcards.strat_values1A}' '{wildcards.strat_values1B}' '{wildcards.strat_name2}' '{wildcards.strat_values2A}' '{wildcards.strat_values2B}' {output}
    """

rule DE_heatmaps:
  input:
    maybe_skip_annotation(do_annot),
    'data/gtf.RDS',
    [results_dir + "DE/{test_name}/{strat_name1}/{test_name}-{value1}-{value2}_DE_Strat1-{strat_name1}-{strat_values1A}-{strat_values1B}_Strat2-{strat_name2}-{strat_values2A}-{strat_values2B}.tsv".format(test_name = comp_dict[dkey]['test_name'], value1 = comp_dict[dkey]['value1'], value2 = comp_dict[dkey]['value2'], strat_name1 = comp_dict[dkey]['strat_name1'], strat_values1A = comp_dict[dkey]['strat_values1A'], strat_values1B = comp_dict[dkey]['strat_values1B'], strat_name2 = comp_dict[dkey]['strat_name2'], strat_values2A = comp_dict[dkey]['strat_values2A'], strat_values2B = comp_dict[dkey]['strat_values2B']) for dkey in comp_dict.keys()]
  params:
    scripts = config['SC_scripts']
  output:
    results_dir + "DE/DE_Heatmap_within_strat.pdf",
    results_dir + "DE/DE_Heatmap_within_test.pdf"
  shell:
    "Rscript {params.scripts}/DE_heatmaps.R {output} {input}"

rule fgsea_test_DE:
  input:
    maybe_skip_annotation(do_annot),
    #results_dir + 'DE/{test}/{strat}/DE_results.tsv',
    results_dir + 'DE/{test_name}/{strat_name1}/{test_name}-{value1}-{value2}_DE_Strat1-{strat_name1}-{strat_values1A}-{strat_values1B}_Strat2-{strat_name2}-{strat_values2A}-{strat_values2B}.tsv',
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
    tsv = results_dir + 'DE/{test_name}/{strat_name1}/{test_name}-{value1}-{value2}_GSEA_Strat1-{strat_name1}-{strat_values1A}-{strat_values1B}_Strat2-{strat_name2}-{strat_values2A}-{strat_values2B}.tsv',
    pdf = results_dir + 'DE/{test_name}/{strat_name1}/{test_name}-{value1}-{value2}_GSEA_Strat1-{strat_name1}-{strat_values1A}-{strat_values1B}_Strat2-{strat_name2}-{strat_values2A}-{strat_values2B}.pdf'
  shell:
    "Rscript {params.scripts}/fgsea.R {input} '{params.gmt}' '{params.species}' {output} '{params.custom_sets}'"

rule fgsea_bubble_plots:
  input:
    maybe_skip_annotation(do_annot),
    'data/gtf.RDS',
    [results_dir + "DE/{test_name}/{strat_name1}/{test_name}-{value1}-{value2}_GSEA_Strat1-{strat_name1}-{strat_values1A}-{strat_values1B}_Strat2-{strat_name2}-{strat_values2A}-{strat_values2B}.tsv".format(test_name = comp_dict[dkey]['test_name'], value1 = comp_dict[dkey]['value1'], value2 = comp_dict[dkey]['value2'], strat_name1 = comp_dict[dkey]['strat_name1'], strat_values1A = comp_dict[dkey]['strat_values1A'], strat_values1B = comp_dict[dkey]['strat_values1B'], strat_name2 = comp_dict[dkey]['strat_name2'], strat_values2A = comp_dict[dkey]['strat_values2A'], strat_values2B = comp_dict[dkey]['strat_values2B']) for dkey in comp_dict.keys()]
  params:
    scripts = config['SC_scripts']
  output:
    results_dir + "DE/GSEA_Bubble_within_strat.pdf",
    results_dir + "DE/GSEA_Bubble_within_test.pdf"
  shell:
    "Rscript {params.scripts}/fgsea_bubble_plots.R {output} {input}"

rule Excel_DE:
  input:
    [results_dir + "DE/{test_name}/{strat_name1}/{test_name}-{value1}-{value2}_DE_Strat1-{strat_name1}-{strat_values1A}-{strat_values1B}_Strat2-{strat_name2}-{strat_values2A}-{strat_values2B}.tsv".format(test_name = comp_dict[dkey]['test_name'], value1 = comp_dict[dkey]['value1'], value2 = comp_dict[dkey]['value2'], strat_name1 = comp_dict[dkey]['strat_name1'], strat_values1A = comp_dict[dkey]['strat_values1A'], strat_values1B = comp_dict[dkey]['strat_values1B'], strat_name2 = comp_dict[dkey]['strat_name2'], strat_values2A = comp_dict[dkey]['strat_values2A'], strat_values2B = comp_dict[dkey]['strat_values2B']) for dkey in comp_dict.keys()]
  params:
    scripts = config['SC_scripts']
  output:
    results_dir + "DE/DE.xls",
    results_dir + "DE/DE_sig.xls"
  shell:
    "Rscript {params.scripts}/Combine_results.R {output} 0.05 True {input}"

rule Excel_fgsea:
  input:
    [results_dir + "DE/{test_name}/{strat_name1}/{test_name}-{value1}-{value2}_GSEA_Strat1-{strat_name1}-{strat_values1A}-{strat_values1B}_Strat2-{strat_name2}-{strat_values2A}-{strat_values2B}.tsv".format(test_name = comp_dict[dkey]['test_name'], value1 = comp_dict[dkey]['value1'], value2 = comp_dict[dkey]['value2'], strat_name1 = comp_dict[dkey]['strat_name1'], strat_values1A = comp_dict[dkey]['strat_values1A'], strat_values1B = comp_dict[dkey]['strat_values1B'], strat_name2 = comp_dict[dkey]['strat_name2'], strat_values2A = comp_dict[dkey]['strat_values2A'], strat_values2B = comp_dict[dkey]['strat_values2B']) for dkey in comp_dict.keys()]
  params:
    scripts = config['SC_scripts']
  output:
    results_dir + "DE/GSEA.xls",
    results_dir + "DE/GSEA_sig.xls"
  shell:
    "Rscript {params.scripts}/Combine_results.R {output} 0.05 True {input}"

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

