import argparse
import logging
import os
import sys
import glob
import itertools
import time

import requests
import pandas as pd
from bs4 import BeautifulSoup as bs

# Set up log
log = logging.getLogger()
formatter = logging.Formatter(
    fmt = '%(asctime)s %(levelname)s %(message)s',
    datefmt = '%m/%d/%Y %I:%M:%S %p')
log.setLevel(logging.DEBUG)

# Log stream handler
sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)
log.addHandler(sh)

# Set up argument parser
desc = ( "Automatically submit jobs to the JGI MiSI tool. Query and reference" +
         " fasta files should contain the predicted genes (DNA) of specified" +
         " genomes. No tRNA or rRNA features should be included. No fasta " +
         " files greater than 500MB in size will be accepted." )
ap = argparse.ArgumentParser(description = desc)
ap.add_argument(
    "qry_dir_path",
    help = "The full path to a directory containing query genomes.")
ap.add_argument(
    "ref_dir_path",
    help = "The full path to a directory containing reference genomes.")
ap.add_argument(
    "out_dir_path",
    help = "The full path to a directory for the output.")
ap.add_argument(
    '--start',
    help = ( "Start submitting jobs at this job number." +
             " Useful for restarting this script." +
             " Defaults to 1." ),
    default = 1,
    metavar = 'NUMBER')
args = ap.parse_args()

# Command line path strings to paths
try:
    assert os.path.normpath(args.qry_dir_path)
    assert os.path.normpath(args.ref_dir_path)
    assert os.path.normpath(args.out_dir_path)
except AssertionError:
    log.error("Non-existant path(s) entered; exiting.")
    sys.exit(1)

# Log file handler
fh = logging.FileHandler(os.path.join(args.out_dir_path, "misi-robot_log.txt"))
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
log.addHandler(fh)

# Accepted fasta extensions
fasta_exts = ('.fasta', '.fna', '.faa', '.fa')

# Set up list of combinations of query and reference fastas
log.debug("Getting list of fasta files from query directory.")
qry_list = []
for qry_file in glob.glob(os.path.join(args.qry_dir_path, "*")):
    has_fasta_ext = qry_file.endswith(fasta_exts)
    file_size = os.path.getsize(qry_file)
    if   has_fasta_ext and file_size <  500000000:
        qry_list.append(qry_file)
    elif has_fasta_ext and file_size >= 500000000:
        log.error("{} exceeds the 500MB submission limit.".format(qry_file) + 
                  " Please remove this file and try again. Exiting.")
        sys.exit(1)
    else:
        log.debug("{} is not a valid fasta file, skipping.".format(qry_file))
        continue

log.debug("Getting list of fasta files from reference directory.")
ref_list = []
for ref_file in glob.glob(os.path.join(args.ref_dir_path, "*")):
    has_fasta_ext = ref_file.endswith(fasta_exts)
    file_size = os.path.getsize(ref_file)
    if   has_fasta_ext and file_size <  500000000:
        ref_list.append(ref_file)
    elif has_fasta_ext and file_size >= 500000000:
        log.error("{} exceeds the 500MB submission limit.".format(ref_file) + 
                  " Please remove this file and try again. Exiting.")
        sys.exit(1)
    else:
        log.debug("{} is not a valid fasta file, skipping.".format(ref_file))
        continue

cmb_list = list(itertools.product(qry_list,ref_list))

tot_subm = len(cmb_list)
log.info("Total jobs to submit: {}".format(tot_subm))

url = "https://ani.jgi-psf.org/html/calc_results.php"
out_file_path = os.path.join(args.out_dir_path, "misi-results.csv")

# Begin submitting data
results = []
for i, cmb in enumerate(cmb_list):
    if i+1 >= args.start:

        # Create dicts for the form
        data  = {"op"    : "seqBasedCalc"}
        files = {"file1" : open(cmb[0],"rb"),
                 "file2" : open(cmb[1],"rb")}

        log.debug("Attempting submission: {}/{} {}".format(i+1, tot_subm, cmb))
        
        # Post the data
        try:
            post = requests.post(url, data = data, files = files)
            rtrn = post.content.decode()
            html = bs(rtrn, 'html.parser')
        except requests.exceptions.RequestException as e:
            log.error("Submission failed: {}/{}".format(i+1, tot_subm))
            log.error("There was an error submitting to the JGI server:")
            log.error(e)
            log.error("You can attempt restarting at this submission. Exiting.")
            raise

        # Convert html tables to pandas tables
        html_tables = html.find_all("table", id="smallTable")
        df_primary = pd.read_html(html_tables[0].prettify())[0]
        df_details = pd.read_html(html_tables[1].prettify(), index_col = 0)[0].T

        # Validate results
        try:
            assert not df_primary.isnull().all().all()
            assert not df_details.isnull().all().all()
        except AssertionError:
            if os.path.basename(cmb[0]) == os.path.basename(cmb[1]): 
                log.info("Submission completed: {}/{}".format(i+1, tot_subm))
                continue # JGI outputs NaNs for identical query-reference pairs
            else:
                log.error("Submission failed: {}/{}".format(i+1, tot_subm))
                log.error("The JGI server unexpectedly returned no results.")
                log.error("You can attempt restarting at this submission. Exiting")
                raise

        # Merge on temp columns
        df_primary['tmp'] = 1
        df_details['tmp'] = 1
        df = pd.merge(df_primary, df_details, on = ['tmp'])
        df = df.drop('tmp', axis = 1)

        # Append to results
        if os.path.exists(out_file_path): 
            df.to_csv(out_file_path, index = False, mode = "a", header=False)
        else:
            df.to_csv(out_file_path, index = False)

        log.info("Submission completed: {}/{}".format(i+1, tot_subm))
        time.sleep(5)

    else:
        log.debug("Skipping {}".format(cmb))
        continue

