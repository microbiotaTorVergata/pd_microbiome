# work dir
PATH_WORK="/../../work_dir"
# path to greengenes
PATH_DB_SEQ="/../../"
PATH_DB_TAX="/../../"
# path to bin (usearch61)
PATH_BIN="/../../bin"
export PATH="$PATH_BIN:$PATH"                      		        
#!/bin/bash
# temp dir
mkdir -p tmpdir 

source activate qiime1

# temp dir
mkdir -p tmpdir 

# qiime config
echo "# qiime_config">>qiime_config
echo "# WARNING: DO NOT EDIT OR DELETE Qiime/qiime/support_files/qiime_config">>qiime_config
echo "# To overwrite defaults, copy this file to $HOME/.qiime_config or a full path">>qiime_config
echo "# specified by $QIIME_CONFIG_FP and edit that copy of the file.">>qiime_config
echo "# Find details on this file at http://qiime.org/install/qiime_config.html">>qiime_config
echo "">>qiime_config
echo "cluster_jobs_fp">>qiime_config
echo "blastmat_dir">>qiime_config
echo "blastall_fp">>qiime_config	
echo "pynast_template_alignment_fp">>qiime_config
echo "pynast_template_alignment_blastdb">>qiime_config
echo "jobs_to_start	1">>qiime_config
echo "seconds_to_sleep	1">>qiime_config
echo "temp_dir	$PATH_WORK/tmpdir">>qiime_config
echo "denoiser_min_per_core	50">>qiime_config
echo "topiaryexplorer_project_dir">>qiime_config
echo "torque_queue	None">>qiime_config
echo "pick_otus_reference_seqs_fp">>qiime_config
echo "assign_taxonomy_reference_seqs_fp">>qiime_config
echo "assign_taxonomy_id_to_taxonomy_fp">>qiime_config
echo "sc_queue	false">>qiime_config
echo "slurm_queue">>qiime_config
echo "slurm_memory">>qiime_config

# validate mapping file
validate_mapping_file.py -m newmapping.txt -o mapping -b -p -s 

source deactivate
