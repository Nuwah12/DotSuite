""" This script will transform given BED files containing features of interest into a 'feature matrix'
These will represent the genome's architectural and functional features in a 2D matrix form

Noah Burget
3/14/23
"""

import cooler
import bioframe
import pybedtools
import argparse
import logging
import pandas as pd
from progress.bar import *

def make_bins(binsize, chromsizes):
    """Make bins of size binsize
    Using bioframe's binnify command
    """ 
    binned_genome = bioframe.binnify(chromsizes = chromsizes, binsize = binsize, rel_ids = True)
    return binned_genome

def main():
    #### Arguments
    parser = argparse.ArgumentParser(description = "Make a feature matrix from a given set of BED files")
    parser.add_argument("--features", "-f", type=str, nargs = "*", help="Path to bed files containing feature regions. Format should be chr,start,end", required=True)
    parser.add_argument("--featureNames", "-n", type=str, nargs = "*", help="Names of columns to be associated with features", required=True)
    parser.add_argument("--resolution", "-r", type=int, help="Size (in bp) of bins. Typically in the range of 400 - 25000", required=True)
    parser.add_argument("--genome", "-g", type=str, help="Genome assembly used to align feature samples", required=False, default='hg38')
    parser.add_argument("--prefix", "-p", type=str, help="Prefix for output file: {prefix}_{res}bp.matrix", default="featureMatrix")
    #parser.add_argument("--dots", "-d", type=str, help="Path to dots (chr1,start1,end1,chr2,start2,end2), whose Anchors will be present as features in the output matrix")
    args=parser.parse_args()

    ### Progress bar
    bar = IncrementalBar('Processing Features', max = len(args.features))

    ### Compatibility checks
    # Feature file array size == feature name array size
    if not len(args.features) == len(args.featureNames):
        raise ValueError("Each feature file must be associated with a name!")

    ### Read in each featureFile
    logging.info("Reading and organizing supplied feature files...")
    featureFiles=[]
    for i in args.features:
        f = bioframe.read_table(i, schema = 'bed3')
        featureFiles.append(f)

    ### Make dict of {featureName : featureFile}
    feature_dict={}
    for i in range(0, len(args.features)):
        feature_dict[args.featureNames[i]] = featureFiles[i]
    
    ### Overlap each featureFile with the binned genome
    # get chromsizes
    chromsizes = bioframe.fetch_chromsizes(args.genome)
    # bin the genome:
    binned_genome = make_bins(int(args.resolution), chromsizes) 
    ### Give bins unique ids 
    binned_genome["rel_id"] = range(0,binned_genome.shape[0])

    ### Pybedtools for intersection
    binned_genome = pybedtools.BedTool.from_dataframe(binned_genome)
    binned_features_master = {}
    for feature in feature_dict:
        binned_features_dict = {}
        feature_bedtool = pybedtools.BedTool.from_dataframe(feature_dict[feature])
        # intersect with binned genome
        binned_features = binned_genome.intersect(feature_bedtool, wao=True)
        # binned_features.to_dataframe(disable_auto_names=True, header=None).to_csv('overlap_test.tsv', sep="\t")
        ### Make column for this feature
        # dict structure: {bin_id : feature_status (1/0)}
        binned_features = binned_features.to_dataframe(disable_auto_names=True, header=None)
        # Loop through entire overlapped dataframe, if overlap size == 0 that bin does not have this feature, else it does 
        # Progress bar
        fulllen_bar = ChargingBar('Binning feature {}'.format(feature), max = len(binned_features.index))
        for i in range(0, len(binned_features.index)):
            if binned_features.iloc[i,7] != 0:
                binned_features_dict[binned_features.iloc[i,3]] = 1
            else:
                binned_features_dict[binned_features.iloc[i,3]] = 0
            fulllen_bar.next()
        binned_features_master[feature] = binned_features_dict
        fulllen_bar.finish()
        bar.next()
        
    bar.finish()

    final = pd.DataFrame.from_dict(binned_features_master)
    final.to_csv('{}_{}bp.matrix'.format(args.prefix, args.resolution), sep = "\t")
    
if __name__ == "__main__":
    main()


