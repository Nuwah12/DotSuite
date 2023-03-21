import yaml
import subprocess
import os
import argparse
from multiprocessing import Process

import cooltools
import cooler
import bioframe

import pandas as pd
import numpy as np

def main():
    # Setup arguments
    parser = argparse.ArgumentParser(prog = "DotCaller",
                                     description = "Call dots using 3 locally-based methods - Dotfinder (cooltools), Mustache, and Chromosight.")
    args = parse_args(parser)
    
    # Argument compatibility checks
    if args.res is None and ".mcool" in args.cool:
        raise SyntaxError("Pass in resolution via --res when using a .mcool file as input.")
    for i in args.tools.split(","):
        if i != "chromosight" and i != "mustache" and i != "dotfinder":
            raise SyntaxError("Invalid tool name passed - check --tools parameter.")
        else:
            continue

    cool = read_cooler(args.cool, args.res)
    config = parse_config(args.config)

    #call tools in parallel
    tools = args.tools.split(",")
    proc = []
    for t in tools:
        if t == "chromosight":
            print("Starting chromosight...")
            p = Process(target = call_chromosight, args = (args.cool, config))
            p.start()
            proc.append(p)
        if t == "mustache":
            print("Starting mustache...")
            p = Process(target = call_mustache, args = (args.cool, config))
            p.start()
            proc.append(p)
        if t == "dotfinder":
            print("Starting Dotfinder...")
            p = Process(target = call_dotfinder, args = (cool, config))
            p.start()
            proc.append(p)
    for p in proc:
        p.join()
            
#def run_in_parallel(*funcs, ):
#    proc = []

#    for f in funcs:
#        p = Process(target = f)
#        p.start()
#        proc.append(p)
    
#    for p in proc:
#        p.join()

def parse_args(parser):
    parser.add_argument('cool', type=str, help='Path to .cool/.mcool file to call dots from')
    parser.add_argument('config', type=str, help='Path to config.yaml file')
    parser.add_argument('--tools', '-t', type=str, help='Comma-separated list of tools to use in dot calling, out of chromosight, mustache, and dotfinder (i.e. --tools chromosight,dotfinder,mustache). Default=chromosight,dotfinder,mustache', default='chromosight,dotfinder,mustache')
    parser.add_argument('--res', '-r', type=int, help='Resolution of the matrix upon which dots will be called. Only required if input cool file is .mcool format. Default=10000')

    return parser.parse_args()

def parse_config(config):
    with open(config, 'r') as config_file:
        f = yaml.safe_load(config_file)
    
    return f

def read_cooler(cool, res):
    # Compatibility check for cooler path
    if not os.path.exists(cool):
        raise OSError("Path to cool file is invalid. Check the path you provided.")
    if ".mcool" in cool:
        cool = cooler.Cooler("{}::resolutions//{}".format(cool, res))
    elif ".cool" in cool:
        cool = cooler.Cooler(cool)

    return cool

def call_dotfinder(cool, config):
    # define all parameters
    expected_value_col = config["df_expected_value_column"]
    weight_name = config["df_weight_name"]
    max_dist = config["df_max_distance"]
    max_nan = config["df_max_nan"]
    num_l_bins = config["df_num_lambda_bins"]
    l_bin_fdr = config["df_lambda_bin_fdr"]
    clust_radius = config["df_clustering_radius"]
    tile_size = config["df_tile_size"]
    threads = config["df_threads"]
    out = config["df_out_file"]
    genome = config["genome"]

    if clust_radius == 0:
        clust_radius = None

    # chromosome sizes
    chromsizes = bioframe.fetch_chromsizes(genome)
    cens = bioframe.fetch_centromeres(genome)
    arms = bioframe.make_chromarms(chromsizes, cens)

    # Subset chromsizes to include only chromosomes in the cooler
    arms = arms.set_index("chrom").loc[cool.chromnames].reset_index()

    # Expected interactions
    expected = cooltools.expected_cis(cool, view_df = arms, nproc = threads)

    # Call dots
    dots = cooltools.dots(
            cool, 
            expected=expected,
            expected_value_col=expected_value_col,
            clr_weight_name=weight_name,
            view_df=arms,
            kernels=None,
            max_loci_separation=max_dist,
            max_nans_tolerated=max_nan,
            n_lambda_bins=num_l_bins,
            lambda_bin_fdr=l_bin_fdr,
            clustering_radius=clust_radius,
            cluster_filtering=None,
            tile_size=tile_size,
            nproc=threads)

    dots.to_csv(out, sep = "\t")

def call_mustache(cool_path, config):
    # define all parameters
    p_threshold = config["m_pThreshold"]
    sparsity_threshold = config["m_sparsityThreshold"]
    resolution = config["resolution"]
    threads = config["m_threads"]
    out = config["m_out_file"]
    path = config["mustache_software"]

    subprocess.run(args = '{} -f {} -r {} -o {} -p {} -pt {} -st {}'.format(path, cool_path, resolution, out, threads, p_threshold, sparsity_threshold), shell=True)

def call_chromosight(cool_path, config):
    # define all parameters
    min_dist = config["cs_min_distance"] 
    max_dist = config["cs_max_distance"]
    perc_zero = config["cs_percentage_zero"]
    threads = config["cs_threads"]
    prefix = config["cs_out_prefix"]
    
    subprocess.run(args = 'chromosight detect --threads {} --min-dist {} --max-dist {} --perc-zero {} {} {}'.format(threads, min_dist, max_dist, perc_zero, cool_path, prefix), shell=True)

if __name__ == "__main__":
    main()    
