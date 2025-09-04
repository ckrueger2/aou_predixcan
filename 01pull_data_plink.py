#!/usr/bin/python

#import libraries
try:
    import pandas as pd
except ImportError:
    import subprocess, sys
    subprocess.run([sys.executable, "-m", "pip", "install", "pandas"], check=True)
    import pandas as pd
    
import os
import argparse
import sys
import hail as hl
import subprocess

#initialize hail
hl.init(default_reference='GRCh38', idempotent=True)

#define bucket to save to
bucket = os.getenv('WORKSPACE_BUCKET')
bucket # gs://fc-secure-bb61452f-d5e2-4d26-9227-6a9444241af8/
out_path = f"{bucket}/data/plink"

mt_wgs_path = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/splitMT/hail.mt"
mt = hl.read_matrix_table(mt_wgs_path)
intervals = ['chr1:108309567-110397919']
filt_mt = hl.filter_intervals(mt, [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in intervals])
filt_mt.show()
hl.export_plink(filt_mt, out_path)
  
