#!/usr/bin/env python3
"""Generate mCH correlation with other genes.

This script will ask for an input folder, which contains files of mCH data for each gene.
Correlation will be computed for other genes and a file named after the input gene ID will be output into the specified
    output folder.
"""
import os
import time
import sys
import timeit
from glob import glob
from multiprocessing.dummy import Pool as ThreadPool

import pandas as pd

df_current = None


def find_target_files(directory):
    """Find a list of mCH data files to process in a directory.

    For now, it simply searches for files ending in '_mCH.txt'.

    Arguments:
        directory (str): Target directory.

    Returns:
        list: str of file names.
        """
    return glob(os.path.join(directory, '*_mCH.txt'))


def process_gene_samples(csv_path):
    """Read into DataFrame and pivot mCH data files.

    Arguments:
        csv_path (str): Target path of gene mCH CSV file.

    Returns:
        pandas.DataFrame: pivoted dataframe.
    """
    # Progress indicator.
    print('.', end='')
    sys.stdout.flush()

    df_gene = pd.read_csv(
        csv_path,
        sep='\t',
        header=None,
        names=['geneId', 'samp', 'original', 'normalized'])
    df_pivot = df_gene.pivot(index='geneId', values='original', columns='samp')
    df_pivot.index.name = None
    return df_pivot


def wizard():
    print(
        'This script generates a correlation matrix for the mCH data of a species.'
    )
    print('Use absolute paths with ending slash please.')
    print('-' * 10)
    print('Input path containing mCH files for that species:')
    input_path = input('> ')
    targets = find_target_files(input_path)
    print('-' * 10 +
          '\nFound {} files. Please ENTER to continue, Ctrl+C to quit.'.format(
              len(targets)))
    input('> ')
    if len(targets) == 0:
        print('Nothing to do.')
        exit()

    print('Gathering sample data', end='')
    pivot_pool = ThreadPool(32)
    pivot_results = pivot_pool.map(process_gene_samples, targets)
    pivot_pool.close()
    pivot_pool.join()

    print('\nConcatenating dataframes...', end='')
    sys.stdout.flush()
    start = timeit.default_timer()
    concat_result = pd.concat(pivot_results)
    stop = timeit.default_timer()
    print('done in {} seconds.'.format(stop - start))
    output_path = 'concat_{}.csv'.format(int(time.time()))
    print('Writing concatenated DFs to {}...'.format(output_path))
    concat_result.to_csv(output_path, sep='\t', na_rep='NA')

    del pivot_results

    print('Calculating correlations...', end='')
    sys.stdout.flush()
    start = timeit.default_timer()
    corr = concat_result.T.corr(method='spearman')
    stop = timeit.default_timer()
    print('done in {} seconds.'.format(stop - start))

    output_path = 'correlations_{}.csv'.format(int(time.time()))
    print('Writing output to {}'.format(output_path))
    concat_result.index.name = None
    corr.to_csv(output_path, sep='\t')


if __name__ == '__main__':
    wizard()
