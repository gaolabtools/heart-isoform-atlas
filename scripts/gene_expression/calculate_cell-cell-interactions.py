import glob
import os
import sys

if len(sys.argv) != 4:
    print("Usage: cpdb_run.py path/to/meta path/to/counts path/for/output")
    sys.exit()

cpdb_file_path = "cpdb/v5.0.0/cellphonedb.zip"
meta_file_path = sys.argv[1]
counts_file_path = sys.argv[2]
out_path = sys.argv[3]

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,
    meta_file_path = meta_file_path,
    counts_file_path = counts_file_path,
    counts_data = 'gene_name',
    output_path = out_path,
    separator = '|',
    threshold = 0.1,
    result_precision = 3,
    debug = False,
    output_suffix = None
)
