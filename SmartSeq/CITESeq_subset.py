
### QC step2, sample filter e.g. entire batches gone wrong
### QC step3, optional imputation of RNA and/or protein values
### QC step4, remove cells based on RNA content
rule RNA_cell_filters:
  input:
    seurat = "data/All_data.RDS",
    samples = sample_sheet,
    qc = QC_specs['step4']
  params:
    scripts = config['scripts_dir'],
    batch_name = config['batch'],
    gtf = config['gtf'],
    results_dir = results_dir
  output:
    results_dir + "{batch}/RNA_cell_filtered.RDS",
    results_dir + "{batch}/RNA_filters.pdf",
    results_dir + "{batch}/Excluded_genes.txt",
    results_dir + "{batch}/Filters_Ncells.txt"
  resources:
    mem_mb=get_mem_mb
  shell:
    """
    mkdir -p {params.results_dir}/{wildcards.batch}/
    Rscript {params.scripts}/RNA_cell_filters.R {input} '{wildcards.batch}' '{params.batch_name}' '{params.gtf}' {output}
    """
### QC step5, optional integration of RNA data

rule combine_QC_outs:
  input:
    [os.path.join(results_dir, Batch, 'Filters_Ncells.txt') for Batch in batchs],
    [os.path.join(results_dir, Batch, 'Excluded_genes.txt') for Batch in batchs]
  params:
    scripts = config['scripts_dir'],
  output:
    results_dir + 'Filters_Ncells.txt',
    results_dir + 'Excluded_genes.txt'
  shell:
    "Rscript {params.scripts}/combine_QC_outs.R {output} {input}"

rule cluster:
  input:
    results_dir + "DSB_normalized_data.RDS",
    results_dir + "Excluded_genes.txt"
  params:
    scripts = config['scripts_dir']
  output:
    temp(results_dir + "Clustered.RDS"),
    results_dir + "Clusters.pdf"
  resources: mem_mb=240000
  shell:
    "Rscript {params.scripts}/Clustering_CITESeq.R {input} {output}"

### QC step 7
checkpoint cluster_filter:
  input:
    results_dir + "Clustered.RDS",
    QC_specs['step7']
  params:
    scripts = config['scripts_dir']
  resources: mem_mb=240000
  output:
    results_dir + "Filtered.RDS",
    results_dir + "Filters.txt",
    results_dir + "Cluster_filters.pdf"
  shell:
    "Rscript {params.scripts}/Cluster_filters_CITESeq.R {input} {output}"

rule post_filter_cluster:
  input:
    results_dir + "Filtered.RDS",
    results_dir + "Excluded_genes.txt"
  params:
    scripts = config['scripts_dir'],
  resources: mem_mb=240000
  output:
    temp(results_dir + "Filtered_clustered.RDS"),
    results_dir + "ADT_clusters_filtered.pdf"
  shell:
    "Rscript {params.scripts}/Clustering_CITESeq.R {input} {output}"

#rule UMAPs_unfiltered:
#  input:
#    "results/Clustered.RDS"
#  params:
#    script = config['scripts_dir'],
#  output:
#    "results/UMAPs_unfiltered.pdf"
#  shell:
#    "Rscript {params.script}/UMAPs.R {input} {output}"

#rule UMAPs:
#  input:
#    "results/Filtered_clustered.RDS"
#  params:
#    script = config['scripts_dir'],
#  output:
#    "results/UMAPs.pdf"
#  shell:
#    "Rscript {params.script}/UMAPs.R {input} {output}"
    

rule annotation:
  input:
    results_dir + "Processed_data.RDS"
  params:
    scripts = config['scripts_dir']
  output:
    results_dir + "Mapped.RDS",
    results_dir + "Mapping.pdf"
  resources: mem_mb=240000
  shell:
    "Rscript {params.scripts}/Annotation.R {input} {output}"

rule Cluster_FindAllMarkers:
  input:
    #results_dir + "Processed_data.RDS",
    maybe_skip_annotation(do_annot),
    results_dir + "Excluded_genes.txt"
  params:
    scripts = config['scripts_dir'],
    results_dir = results_dir
  output:
    results_dir + "Cluster_DE/One-all/{assay}-{cluster}.RDS"
  resources: mem_mb=240000
  shell:   
    """
    mkdir -p {params.results_dir}/Cluster_DE/One-all/
    Rscript {params.scripts}/Cluster_DE_once.R {input} {output} {wildcards.assay} {wildcards.cluster}
    """

rule Cluster_FindMarkers:
  input:
    maybe_skip_annotation(do_annot),
    results_dir + "Excluded_genes.txt"
  params:
    scripts = config['scripts_dir'],
    results_dir = results_dir
  output:
    results_dir + "Cluster_DE/One-one/{assay}_{clusters}.RDS"
  resources: mem_mb=240000
  shell:
    """
    mkdir -p {params.results_dir}/Cluster_DE/One-one/
    Rscript {params.scripts}/Cluster_FindMarkers.R {input} {output} {wildcards.assay} {wildcards.clusters}
    """

rule Combine_FindMarkers:
  input:
    expand(os.path.join(results_dir, 'Cluster_DE', 'One-one', "{assay}_{clusters}.RDS"), assay = assays_list, clusters = cluster_combos)
  output:
    results_dir + "Cluster_DE/Singles.txt"
  shell:
    "ls {input} > {output}"

rule SubCluster:
  input:
    maybe_skip_annotation(do_annot),
    results_dir + "Excluded_genes.txt"
  params:
    scripts = config['scripts_dir'],
    results_dir = results_dir
  output:
    results_dir + "SubClusters/{cluster}.RDS"
  resources: mem_mb=240000
  shell:
    """
    mkdir -p {params.results_dir}/SubCluster/
    Rscript {params.scripts}/SubCluster.R {input} {wildcards.cluster} {output}
    """

rule Combine_SubClusters:
  input:
    expand(os.path.join(results_dir, 'SubClusters', "{cluster}.RDS"), cluster = sub_clusters)
  output:
    results_dir + "SubClusters/SubClusters.csv"
  shell:
    "ls {input} {output}"

rule App_data:
  input:
    #results_dir + "Processed_data.RDS",
    #maybe_skip_annotation(config['species'], config['tissue']),
    results_dir + "Cluster_DE/Singles.txt",
    maybe_skip_annotation(do_annot),
    expand(os.path.join(results_dir, 'Cluster_DE', 'One-all', "{assay}-{cluster}.RDS"), assay = assays_list, cluster = clusters_list),
    expand(os.path.join(results_dir, 'SubClusters', "{cluster}.RDS"), assay = assays_list, cluster = sub_clusters)
  params:
    scripts = config['scripts_dir'], 
    project = config['project'],
    username = config['username'],
    server = config['server'],
    results_dir = results_dir
  output:
    results_dir + "App/Data.RData"
  shell:
    """
    mkdir -p {params.results_dir}/App/
    cp {params.scripts}/App.R {params.results_dir}/App/
    cp {params.scripts}/../../sc_functions.R {params.results_dir}/App/
    Rscript {params.scripts}/Create_app_data.R {output} '{params.project}' '{params.username}' '{params.server}' {input}
    """

