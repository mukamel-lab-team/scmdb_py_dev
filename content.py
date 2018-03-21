"""Functions used to generate content. """
import datetime
import csv
import glob
import json
import math
import os
import sys
import time

from collections import OrderedDict
from itertools import groupby, chain
from random import sample

import colorsys
import colorlover as cl
from flask import Blueprint, current_app
from sqlalchemy import exc, text
import numpy as np
from numpy import nan, linspace, arange, random
import pandas as pd
import plotly
from plotly import tools
from plotly.graph_objs import Layout, Annotation, Box, Scatter, Scattergl, Scatter3d, Heatmap
import sqlite3
from sqlite3 import Error

from . import cache, db

content = Blueprint('content', __name__) # Flask "bootstrap"

cluster_annotation_order = ['mL2/3', 'mL4', 'mL5-1', 'mL5-2', 'mDL-1', 'mDL-2', \
                            'mL6-1', 'mL6-2', 'mDL-3', 'mVip', 'mNdnf-1', \
                            'mNdnf-2', 'mPv', 'mSst-1', 'mSst-2', 'None']
methylation_types_order = ['mCH', 'mCG', 'mCA', 'mCHmCG', 'mCHmCA', 'mCAmCG']

class FailToGraphException(Exception):
    """Fail to generate data or graph due to an internal error."""
    pass


@content.route('/content/metadata/')
def get_metadata():

    result = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM cells;").fetchall()
    result = [dict(r) for r in result]

    return json.dumps({"data": result})


@content.route('/content/ensemble_list_dt')
def get_ensemble_list_dt():

    ensemble_list = []
    ensemble_list = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles").fetchall()

    ensemble_dict = [{"ensemble_id": d.ensemble_id,
                      "ensemble_name": d.ensemble_name,
                      "total_cells": d.datasets} for d in ensemble_list]

    return_json = json.dumps(ensemble_dict)
    return return_json 

@content.route('/content/ensemble_list')
def get_ensemble_list():
    
    ensemble_list=[]
    ensemble_list = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles").fetchall()

    total_cell_each_dataset = db.get_engine(current_app, 'methylation_data').execute("SELECT dataset, COUNT(*) as `num` FROM cells GROUP BY dataset").fetchall()
    total_cell_each_dataset = [ {d['dataset']: d['num']} for d in total_cell_each_dataset ]
    total_cell_each_dataset = { k.split('_',maxsplit=1)[1]: v for d in total_cell_each_dataset for k, v in d.items() }

    ensembles_cell_counts = []
    for ensemble in ensemble_list:
        ensemble_tbl = 'Ens' + str(ensemble['ensemble_id'])
        query = "SELECT dataset, COUNT(*) as `num` FROM cells INNER JOIN ? ON cells.cell_id = ?.cell_id GROUP BY dataset"
        cell_counts = db.get_engine(current_app, 'methylation_data').execute(query, (ensemble,ensemble,)).fetchall()
        cell_counts = [ {d['dataset']: d['num']} for d in cell_counts]
        cell_counts = { k.split('_',maxsplit=1)[1]: v for d in cell_counts for k, v in d.items() }
        ensembles_cell_counts.append( {"id": ensemble['ensemble_id'], "ensemble": ensemble['ensemble_name'], "ens_counts": cell_counts} )
    ensembles_json_list = []
    for ens in ensembles_cell_counts:
        total_cells = 0
        datasets_in_ensemble = []
        ens_dict = {}
        for dataset, count in ens['ens_counts'].items():
            total_cells = total_cells + count
            datasets_in_ensemble.append(dataset)
            ens_dict[dataset] = str(count) + '/' + str(total_cell_each_dataset[dataset])
        ens_dict["id"] = ens["id"]
        ens_dict["ensemble"] = ens['ensemble']
        ens_dict["datasets"] = "\n".join(datasets_in_ensemble)
        ens_dict["total_cells"] = total_cells

        ens_regions_query = "SELECT DISTINCT(ABA_name) FROM ABA_regions WHERE `code` IN (" + ",".join("'"+x.split('_')[0]+"'" for x in datasets_in_ensemble) + ")"
        ens_regions_result = db.get_engine(current_app, 'methylation_data').execute(ens_regions_query).fetchall()
        ens_regions_result = [d['ABA_name'] for d in ens_regions_result]
        ens_dict["regions"] = ", ".join(ens_regions_result)
        ensembles_json_list.append(ens_dict)

    all_datasets = db.get_engine(current_app, 'methylation_data').execute("SELECT DISTINCT(dataset) FROM cells").fetchall()
    all_datasets = [d['dataset'] for d in all_datasets]
    all_datasets = [d.split('_', maxsplit=1)[1] for d in all_datasets]

    ## Matches code in dataset name with Allen Brain Atlas regions. ex: 'CEMBA_3C_171206' -> '3C' -> 'MOp'(Primary Motor Cortex)
    dataset_regions = [d.split('_', maxsplit=1)[0] for d in all_datasets]
    zipped_regions = list(zip(all_datasets, dataset_regions))
    ABA_regions = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ABA_regions WHERE `code` IN (" + ",".join("'"+x+"'" for x in dataset_regions) + ")").fetchall()
    ABA_regions = [tuple(d) for d in ABA_regions]

    datasets_with_ABA_list = []
    for dataset, region in zipped_regions:
        for code, ABA_name in ABA_regions:
            if region == code:
                dataset_dict = {"dataset": dataset, "region": ABA_name}
                break
        datasets_with_ABA_list.append(dataset_dict)

    data_dict = {"columns": datasets_with_ABA_list,
                 "data":ensembles_json_list}
    ens_json = json.dumps(data_dict)

    return ens_json


# Utilities
@cache.memoize(timeout=1800)
def ensemble_exists(ensemble):
    """Check if data for a given ensemble exists by looking for its data directory.

    Arguments:
        ensemble (str): Name of ensemble.

    Returns:
        bool: Whether if given ensemble exists
    """

    return db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles WHERE ensemble_name=%s", (ensemble,)).fetchone() != None


@cache.memoize(timeout=1800)
def gene_exists(ensemble, methylation_type, gene):
    """Check if data for a given gene of ensemble exists by looking for its data directory.

    Arguments:
        ensemble (str): Name of ensemble.
        methylation_type (str): Type of methylation to visualize. "mch" or "mcg"
        gene (str): Ensembl ID of gene for that ensemble.

    Returns:
        bool: Whether if given gene exists
    """

    gene_table_name = 'gene_' + gene.replace(".", "_")
    return len(db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM information_schema.tables WHERE table_name = %s", (gene_table_name,)).fetchall()) > 0


def build_hover_text(labels):
    """Build HTML for Plot.ly graph labels.

        Arguments:
            labels (dict): Dictionary of attributes to be displayed.

        Returns:
            str: Generated HTML for labels.

        Example:
            >>> build_hover_text({'Test1': 'Value1', 'Example2': 'Words2'})
            'Test1: Value1<br>Example2: Words2'

    """
    text = str()
    for k, v in labels.items():
        text += '{k}: {v}<br>'.format(k=k, v=str(v))

    return text.strip('<br>')


def generate_cluster_colors(num, grouping):
    """Generate a list of colors given number needed.

    Arguments:
        num (int): Number of colors needed. n <= 35.

    Returns:
        list: strings containing RGB-style strings e.g. rgb(255,255,255).
    """

    if grouping == 'dataset' and num > 2 and num <= 9:
        c = cl.scales[str(num)]['qual']['Set1'] 
        return c

    if num>18:
        # c = ['hsl('+str(round(h*1.8 % 360))+',50%,50%)' for h in linspace(0, 360, num)]
        c = ['rgb'+str(colorsys.hls_to_rgb((h*1.8/360), 0.5, 0.5)) for h in linspace(0, 360, num)]
    else:
        # c = ['hsl('+str(round(h*1.3 % 360))+',50%,50%)' for h in linspace(0, 360, num)]
        c = ['rgb'+str(colorsys.hls_to_rgb((h*1.3 / 360), 0.5, 0.5)) for h in linspace(0, 360, num)]

    c=c+c
    return c


# def randomize_cluster_colors(num):
#     """Generates random set of colors for tSNE cluster plot.
# 
#     Arguments:
#         None
# 
#     Returns:
#         list: dict items.
#             'colors' = new set of colors for each trace in rgb.
#             'num_colors' = number of colors to be used
#             'cluster_color_#' = indexes of traces to be assigned the new color
# 
#     """
#     cache.delete_memoized(generate_cluster_colors)
#     try:
#         new_colors = {'colors': generate_cluster_colors(num)}
#         new_colors['num_colors'] = num
#         new_colors.update(trace_colors)
#         return new_colors
#     except NameError:
#         time.sleep(2)
#         randomize_cluster_colors(num)
# 

def set_color_by_percentile(this, start, end):
    """Set color below or above percentiles to their given values.

    Since the Plot.ly library handles coloring, we work directly with mCH values in this function. The two percentiles
    are generated by the pd library from the plot-generating method.

    Arguments:
        this (float): mCH value to be compared.
        start (float): Lower end of percentile.
        end (float): Upper end of percentile.

    Returns:
        int: Value of `this`, if it is within percentile limits. Otherwise return one of two percentiles.
    """
    if str(this) == 'nan':
        return 'grey'
    if this < start:
        return start
    elif this > end:
        return end
    return this


def find_orthologs(mmu_gene_id=str(), hsa_gene_id=str()):
    """Find orthologs of a gene.

    Either hsa_gene_id or mmu_gene_id should be completed.

    Arguments:
        mmu_gene_id (str): Ensembl gene ID of mouse.
        hsa_gene_id (str): Ensembl gene ID of human.

    Returns:
        dict: hsa_gene_id and mmu_gene_id as strings.
    """
    if not mmu_gene_id and not hsa_gene_id:  # Should have at least one.
        return {'mmu_gene_id': None, 'hsa_gene_id': None}

    conn = sqlite3.connect(
        '{}/datasets/orthologs.sqlite3'.format(current_app.config['DATA_DIR']))

    # This ensures dictionaries are returned for fetch results.
    conn.row_factory = sqlite3.Row  

    cursor = conn.cursor()
    query_key = 'mmu_gene_id' if mmu_gene_id else 'hsa_gene_id'
    query_value = mmu_gene_id or hsa_gene_id
    query_value = query_value.split('.')[0] 

    try:
        cursor.execute("SELECT * FROM orthologs WHERE ? LIKE ?", (query_key, '%'+query_value+'%',))
    except sqlite3.Error as e:
        now = datetime.datetime.now()
        print("[{}] ERROR in app(find_orthologs): {}".format(str(now), e))
        sys.stdout.flush()
        return {'mmu_gene_id': None, 'hsa_gene_id': None}

    query_results = cursor.fetchone()
    if not query_results:
        return {'mmu_gene_id': None, 'hsa_gene_id': None}
    else:
        return dict(query_results)


@cache.cached(timeout=3600)
def all_gene_modules():
    """Generate list of gene modules for populating gene modules selector.
    Arguments:
        None
    Returns:
        list of gene module names. 
    """

    modules_result = db.get_engine(current_app, 'methylation_data').execute("SELECT DISTINCT(module) FROM gene_modules").fetchall()
    modules = [{'module': module['module']} for module in modules_result]

    return modules


@cache.memoize(timeout=1800)
def get_genes_of_module(ensemble, module):
    """Generates list of genes in selected module.
    Arguments:
        module (str): Name of module to query for.
    Returns:
        Dataframe of gene_name and gene_id of each gene in the module for the corresponding
        ensemble.
    """

    modules_result = db.get_engine(current_app, 'methylation_data').execute("SELECT module, mmu_gene_id, mmu_gene_name FROM gene_modules WHERE module=%s", (module,)).fetchall()
    genes_in_module = [ {'module': d['module'], 'gene_id': d['mmu_gene_id'], 'gene_name': d['mmu_gene_name']} for d in modules_result ]

    return genes_in_module


@cache.memoize(timeout=3600)
def median_cluster_mch(gene_info, grouping, clustering):
    """Returns median mch level of a gene for each cluster.

        Arguments:
            gene_info (dict): mCH data for each sample. Keys are samp(cell), tsne_x, tsne_y, cluster_label, cluster_ordered, original, normalized.
            level (str): "original" or "normalized" methylation values.

        Returns:
            dict: Cluster_label (key) : median mCH level (value).
    """

    if grouping == 'annotation':
        gene_info.fillna({'annotation_'+clustering: 'None'}, inplace=True)
    return gene_info.groupby(grouping+'_'+clustering, sort=False)[gene_info.columns[-2]].median()


@cache.memoize(timeout=3600)
def snATAC_data_exists(ensemble_id):
    ensemble = 'Ens' + str(ensemble_id)
    result = db.get_engine(current_app, 'snATAC_data').execute("SHOW TABLES LIKE %s", (ensemble,)).fetchone()

    if result is None:
        return 0
    else:
        return 1


@cache.memoize(timeout=3600)
def get_ensemble_info(ensemble_name=str(), ensemble_id=str()):
    """
    Gets information regarding an ensemble. Requires either the ensemble name or ensemble id
    """
    
    if ensemble_name:
        result = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles WHERE ensemble_name=%s", (ensemble_name,)).fetchone()
    else:
        ensemble_id = int(filter(str.isdigit, ensemble_id))
        result = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles WHERE ensemble_id=%s", (ensemble_id,)).fetchone()

    return result


@cache.memoize(timeout=3600)
def get_methylation_tsne_options(ensemble):
    """
    Get all available options for tsne plot for selected ensemble.
    """

    if ";" in ensemble: # Prevent SQL injection since table names aren't parameterizable
        return None

    query = "SELECT * FROM %(ensemble)s LIMIT 1" % {'ensemble': ensemble,}
    
    try:
        df = pd.read_sql(query, db.get_engine(current_app, 'methylation_data'))
    except exc.ProgrammingError as e:
        now = datetime.datetime.now()
        print("[{}] Error in app(get_methylation_tsne_options): {}".format(str(now), e))
        sys.stdout.flush()
        return None

    df_tsne = df.filter(regex='^tsne_x_', axis='columns')
    df_cluster = df.filter(regex='^cluster_', axis='columns')

    list_tsne_types = [x.split('tsne_x_')[1] for x in df_tsne.columns.values]
    list_mc_types_tsne = sorted(list(set([x.split('_')[0] for x in list_tsne_types])), key=lambda mC_type: methylation_types_order.index(mC_type))
    list_dims_tsne_first = sorted(list(set([int(x.split('_')[1].replace('ndim','')) for x in list_tsne_types if list_mc_types_tsne[0] == x.split('_')[0]])))
    list_perp_tsne_first = sorted(list(set([int(x.split('_')[2].replace('perp', '')) for x in list_tsne_types if (list_mc_types_tsne[0]+'_ndim'+str(list_dims_tsne_first[0])) == (x.split('_')[0] +'_'+ x.split('_')[1])])))

    list_clustering_types = [x.split('cluster_')[1] for x in df_cluster.columns.values]

    #generate query for getting number clusters for each clustering type
    num_clusters_query = "SELECT "
    for i, clustering in enumerate(list_clustering_types):
        num_clusters_query += "MAX(cluster_{0}) as {0}, ".format(clustering)
    num_clusters_query = num_clusters_query[:-2] #Gets rid of last ", " which causes a MySQL syntax error
    num_clusters_query += " FROM {}".format(ensemble)

    result = db.get_engine(current_app, 'methylation_data').execute(num_clusters_query).fetchone()

    dict_clustering_types_and_numclusters = OrderedDict()
    for clustering_type, num_clusters in zip(list_clustering_types, result):
        dict_clustering_types_and_numclusters[clustering_type] = num_clusters

    list_mc_types_clustering = sorted(list(set([x.split('_')[0] for x in list_clustering_types])), key=lambda mC_type: methylation_types_order.index(mC_type))
    list_algorithms_clustering = sorted(list(set([x.split('_')[1] for x in list_clustering_types if list_mc_types_clustering[0] == x.split('_')[0]])))
    list_npc_clustering = sorted(list(set([int(x.split('_')[2].replace('npc', '')) for x in list_clustering_types if (list_mc_types_clustering[0]+'_'+list_algorithms_clustering[0]) in x])))
    list_k_clustering = sorted(list(set([int(x.split('_')[3].replace('k', '')) for x in list_clustering_types if (list_mc_types_clustering[0]+'_'+list_algorithms_clustering[0]+'_npc'+str(list_npc_clustering[0])) in x])))

    return {'all_tsne_types': list_tsne_types, 
            'tsne_methylation': list_mc_types_tsne,
            'tsne_dimensions': list_dims_tsne_first,
            'tsne_perplexity': list_perp_tsne_first,
            'all_clustering_types': list_clustering_types,
            'all_clustering_types2': dict_clustering_types_and_numclusters,
            'clustering_methylation': list_mc_types_clustering,
            'clustering_algorithms': list_algorithms_clustering,
            'clustering_npc': list_npc_clustering,
            'clustering_k': list_k_clustering,}


@cache.memoize(timeout=3600)
def get_snATAC_tsne_options(ensemble):
    """
    Get all available options for tsne plot for selected ensemble.
    """

    if ";" in ensemble: # Prevent SQL injection since table names aren't parameterizable
        return None

    query = "SELECT * FROM %(ensemble)s LIMIT 1" % {'ensemble': ensemble,}
    
    try:
        df = pd.read_sql(query, db.get_engine(current_app, 'snATAC_data'))
    except exc.ProgrammingError as e:
        now = datetime.datetime.now()
        print("[{}] Error in app(get_snATAC_tsne_options): {}".format(str(now), e))
        sys.stdout.flush()
        return None

    df_tsne = df.filter(regex='^tsne_x_', axis='columns')
    df_cluster = df.filter(regex='^cluster_', axis='columns')

    list_tsne_types = [x.split('tsne_x_')[1] for x in df_tsne.columns.values]
    list_dims_tsne_first = sorted(list(set([int(x.split('_')[1].replace('ndim','')) for x in list_tsne_types])))
    list_perp_tsne_first = sorted(list(set([int(x.split('_')[2].replace('perp', '')) for x in list_tsne_types])))

    list_clustering_types = [x.split('cluster_')[1] for x in df_cluster.columns.values]

    #generate query for getting number clusters for each clustering type
    num_clusters_query = "SELECT "
    for i, clustering in enumerate(list_clustering_types):
        num_clusters_query += "MAX(cluster_{0}) as {0}, ".format(clustering)
    num_clusters_query = num_clusters_query[:-2] #Gets rid of last ", " which causes a MySQL syntax error
    num_clusters_query += " FROM {}".format(ensemble)

    result = db.get_engine(current_app, 'snATAC_data').execute(num_clusters_query).fetchone()

    dict_clustering_types_and_numclusters = OrderedDict()
    for clustering_type, num_clusters in zip(list_clustering_types, result):
        dict_clustering_types_and_numclusters[clustering_type] = num_clusters

    list_algorithms_clustering = sorted(list(set([x.split('_')[1] for x in list_clustering_types])))
    list_npc_clustering = sorted(list(set([int(x.split('_')[2].replace('npc', '')) for x in list_clustering_types])))
    list_k_clustering = sorted(list(set([int(x.split('_')[3].replace('k', '')) for x in list_clustering_types])))

    return {'all_tsne_types': list_tsne_types, 
            'tsne_dimensions': list_dims_tsne_first,
            'tsne_perplexity': list_perp_tsne_first,
            'all_clustering_types': list_clustering_types,
            'all_clustering_types2': dict_clustering_types_and_numclusters,
            'clustering_algorithms': list_algorithms_clustering,
            'clustering_npc': list_npc_clustering,
            'clustering_k': list_k_clustering,}


@cache.memoize(timeout=3600)
def get_gene_by_name(ensemble, gene_query):
    """Match gene names of a ensemble.

    Arguments:
        ensemble (str): Name of ensemble.
        gene_query (str): Query string of gene name.

    Returns:
        DataFrame: Info for queried gene. Columns are gene_id, gene_name, chr, start, end, strand, gene_type.
    """

    gene_query = gene_query.lower()

    df = pd.read_sql("SELECT * FROM genes WHERE lower(gene_name) LIKE %s", params=(gene_query+"%",), con=db.get_engine(current_app, 'methylation_data'))

    return df.to_dict('records')
    

@cache.memoize(timeout=3600)
def get_gene_by_id(ensemble, gene_query):
    """Match gene ID of a ensemble.

    Arguments:
        ensemble (str): Name of ensemble.
        gene_query (str): Query string of gene ID.

    Returns:
        DataFrame: Info for queried gene. Columns are gene_id, gene_name, chr, start, end, strand, gene_type.
    """

    result = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM genes WHERE gene_id LIKE %s", (gene_query+"%",)).fetchone()

    if result is None:
        return {}
    else:
        return dict(result)


def convert_gene_id_mmu_hsa(ensemble, gene_id):
    """Converts mmu gene_ids to hsa gene_ids and vice versa if the ID's do not correspond to
    the ensemble being viewed. Necessary due to caching of last viewed gene IDs.

    Arguments:
        ensemble(str): Name of ensemble.
        gene_id(str): Ensembl ID of gene.of

    Returns:
        str: correct gene_id that corresponds with the ensemble.
    """
    
    if ensemble.startswith('Ens'):
        if "ENSMUSG" not in gene_id:   #ENSUMG is the Ensembl header for mouse genes
            mmu = find_orthologs(hsa_gene_id=gene_id)['mmu_gene_id']
            return find_orthologs(hsa_gene_id=gene_id)['mmu_gene_id']
        else:
            return gene_id
    else:
        if "ENSMUSG" in gene_id:
            gene2 = find_orthologs(mmu_gene_id=gene_id)['hsa_gene_id']
            return find_orthologs(mmu_gene_id=gene_id)['hsa_gene_id']
        else:
            return gene_id


def get_corr_genes(ensemble,query):
    """Get correlated genes of a certain gene of a ensemble. 
    
        Arguments:
            ensemble(str): Name of ensemble.
            query(str): Query string of gene ID.
        
        Returns:
            dict: information of genes that are correlated with target gene.
    """
    db_location = '{}/ensembles/{}/top_corr_genes.sqlite3'.format(
        current_app.config['DATA_DIR'], ensemble)
    try:
        conn = sqlite3.connect(db_location)
    except sqlite3.Error as e:
        now = datetime.datetime.now()
        print("[{}] ERROR in app(get_corr_genes): {}".format(str(now), e))
        sys.stdout.flush()
        return(1)

    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    try:
        cursor.execute('SELECT Gene2, Correlation FROM corr_genes WHERE Gene1 LIKE ? ORDER BY Correlation DESC LIMIT 50', (query + '%',))
    except sqlite3.Error:
        now = datetime.datetime.now()
        print("[{}] ERROR in app(get_corr_genes): Could not load top_corr_genes.sqlite3 for {}".format(str(now), ensemble))
        sys.stdout.flush()
        return(1)

    query_results = list(cursor.fetchall())
    table_data=[]
    for rank, item in enumerate(query_results, 1):
        gene = dict(item)
        geneInfo = get_gene_by_id(ensemble, gene['Gene2'])
        geneInfo['Rank'] = rank
        geneInfo['Corr'] = gene['Correlation']
        table_data.append(geneInfo)
    return table_data


@cache.memoize(timeout=3600)
def get_gene_methylation(ensemble, methylation_type, gene, grouping, clustering, level, outliers, tsne_type='mCH_ndim2_perp20'):
    """Return mCH data points for a given gene.

    Data from ID-to-Name mapping and tSNE points are combined for plot generation.

    Arguments:
        ensemble (str): Name of ensemble.
        methylation_type (str): Type of methylation to visualize. "mCH", "mCG", or "mCA"
        gene (str): Ensembl ID of gene.
        grouping (str): Variable for grouping cells. "cluster", "annotation", or "dataset".
        clustering (str): Different clustering algorithms and parameters. 'lv' = Louvain clustering.
        level (str): "original" or "normalized" methylation values.
        outliers (bool): Whether if outliers should be kept.
        tsne_type (str): Options for calculating tSNE. ndims = number of dimensions, perp = perplexity.

    Returns:
        DataFrame
    """

    # Prevent SQL injected since column names cannot be parameterized.
    if ";" in ensemble or ";" in methylation_type or ";" in grouping or ";" in clustering or ";" in tsne_type:
        return None

    # This query is just to fix gene id's missing the ensemble version number. 
    # Necessary because the table name must match exactly with whats on the MySQL database.
    # Ex. ENSMUSG00000026787 is fixed to ENSMUSG00000026787.3 -> gene_ENSMUSG00000026787_3 (table name in MySQL)
    result = db.get_engine(current_app, 'methylation_data').execute("SELECT gene_id FROM genes WHERE gene_id LIKE %s", (gene+"%",)).fetchone()
    gene_table_name = 'gene_' + result.gene_id.replace('.','_')

    context = methylation_type[1:]

    if 'ndim2' in tsne_type:
        query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, cells.global_%(methylation_type)s, \
            %(ensemble)s.annotation_%(clustering)s, %(ensemble)s.cluster_%(clustering)s, \
            %(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, \
            %(gene_table_name)s.%(methylation_type)s, %(gene_table_name)s.%(context)s \
            FROM cells \
            INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
            LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble,
                                                                                                     'gene_table_name': gene_table_name,
                                                                                                     'tsne_type': tsne_type,
                                                                                                     'methylation_type': methylation_type,
                                                                                                     'context': context,
                                                                                                     'clustering': clustering,}
    else:
        query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, cells.global_%(methylation_type)s, \
            %(ensemble)s.annotation_%(clustering)s, %(ensemble)s.cluster_%(clustering)s, \
            %(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, %(ensemble)s.tsne_z_%(tsne_type)s, \
            %(gene_table_name)s.%(methylation_type)s, %(gene_table_name)s.%(context)s \
            FROM cells \
            INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
            LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble, 
                                                                                                     'gene_table_name': gene_table_name,
                                                                                                     'tsne_type': tsne_type,
                                                                                                     'methylation_type': methylation_type,
                                                                                                     'context': context,
                                                                                                     'clustering': clustering,}
    try:
        df = pd.read_sql(query, db.get_engine(current_app, 'methylation_data'))
    except exc.ProgrammingError as e:
        now = datetime.datetime.now()
        print("[{}] ERROR in app(get_gene_methylation): {}".format(str(now), e))
        sys.stdout.flush()
        return None
    
    if df[context].isnull().all(): # If no data in column, return None 
        return None

    if level == 'original':
        df[methylation_type + '/' + context + '_' + level] = df[methylation_type] / df[context]
    else:
        df[methylation_type + '/' + context + '_' + level] = (df[methylation_type] / df[context]) / df['global_'+methylation_type]

    if not outliers:
        # Outliers not wanted, remove rows > 99%ile
        three_std_dev = df[methylation_type + '/' + context + '_' + level].quantile(0.99) 
        df = df[df[methylation_type + '/' + context + '_' + level] < three_std_dev] 

    
    if grouping == 'annotation':
        df.fillna({'annotation_'+clustering: 'None'}, inplace=True)
        df['annotation_cat'] = pd.Categorical(df['annotation_'+clustering], cluster_annotation_order)
        return df.sort_values(by='annotation_cat')
    elif grouping == 'cluster':
        return df.sort_values(by='cluster_'+clustering)
    else:
        return df


def get_mult_gene_methylation(ensemble, methylation_type, genes, grouping, clustering, level, tsne_type='mCH_ndim2_perp20'):
    """Return averaged methylation data ponts for a set of genes.

    Data from ID-to-Name mapping and tSNE points are combined for plot generation.

    Arguments:
        ensemble (str): Name of ensemble.
        methylation_type (str): Type of methylation to visualize. "mch" or "mcg"
        genes ([str]): List of gene IDs to query.
        grouping (str): Variable for grouping cells. "cluster", "annotation", or "dataset".
        clustering (str): Different clustering algorithms and parameters. 'lv' = Louvain clustering.
        level (str): "original" or "normalized" methylation values.
        tsne_type (str): Options for calculating tSNE. ndims = number of dimensions, perp = perplexity.

    Returns:
        DataFrame
    """

    # Prevent SQL injected since column names cannot be parameterized.
    if ";" in ensemble or ";" in methylation_type or ";" in grouping or ";" in clustering or ";" in tsne_type:
        return None

    context = methylation_type[1:]
    genes = [gene+"%" for gene in genes]

    # This query is just to fix gene id's missing the ensemble version number. 
    # Necessary because the table name must match exactly with whats on the MySQL database.
    # Ex. ENSMUSG00000026787 is fixed to ENSMUSG00000026787.3
    first_query = "SELECT gene_id FROM genes WHERE gene_id LIKE %s" + " OR gene_id LIKE %s" * (len(genes)-1)
    result = db.get_engine(current_app, 'methylation_data').execute(first_query, (genes,)).fetchall()

    gene_table_names = ['gene_' + gene_id[0].replace('.','_') for gene_id in result]

    df_all = pd.DataFrame()
    
    first = True
    for gene_table_name in gene_table_names:
        if 'ndim2' in tsne_type:
            query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, cells.global_%(methylation_type)s, \
                %(ensemble)s.annotation_%(clustering)s, %(ensemble)s.cluster_%(clustering)s, \
                %(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, \
                %(gene_table_name)s.%(methylation_type)s, %(gene_table_name)s.%(context)s \
                FROM cells \
                INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
                LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble, 
                                                                                                         'gene_table_name': gene_table_name,
                                                                                                         'tsne_type': tsne_type,
                                                                                                         'methylation_type': methylation_type,
                                                                                                         'context': context,
                                                                                                         'clustering': clustering,}
        else:
            query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, cells.global_%(methylation_type)s, \
                %(ensemble)s.annotation_%(clustering)s, %(ensemble)s.cluster_%(clustering)s, \
                %(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, %(ensemble)s.tsne_z_%(tsne_type)s, \
                %(gene_table_name)s.%(methylation_type)s, %(gene_table_name)s.%(context)s \
                FROM cells \
                INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
                LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble, 
                                                                                                         'gene_table_name': gene_table_name,
                                                                                                         'tsne_type': tsne_type,
                                                                                                         'methylation_type': methylation_type,
                                                                                                         'context': context,
                                                                                                         'clustering': clustering,}
        try:
            df_all = df_all.append(pd.read_sql(query, db.get_engine(current_app, 'methylation_data')))
        except exc.ProgrammingError as e:
            now = datetime.datetime.now()
            print("[{}] ERROR in app(get_mult_gene_methylation): {}".format(str(now), e))
            sys.stdout.flush()
            return None
        
        if first:
            df_coords = df_all
        first = False

    df_avg_methylation = df_all.groupby(by='cell_id', as_index=False)[[methylation_type, context]].mean()
    df_coords.update(df_avg_methylation)

    if df_coords[context].isnull().all(): # If no data in column, return None 
        return None
    else:
        if level == 'original':
            df_coords[methylation_type + '/' + context + '_' + level] = df_coords[methylation_type] / df_coords[context]
        else:
            df_coords[methylation_type + '/' + context + '_' + level] = (df_coords[methylation_type] / df_coords[context]) / df_coords['global_'+methylation_type]

    if grouping == 'annotation':
        df_coords.fillna({'annotation_'+clustering: 'None'}, inplace=True)
        df_coords['annotation_cat'] = pd.Categorical(df_coords['annotation_'+clustering], cluster_annotation_order)
        return df_coords.sort_values(by='annotation_cat')
    elif grouping == 'cluster':
        return df_coords.sort_values(by='cluster_'+clustering)
    else:
        return df_coords


@cache.memoize(timeout=3600)
def get_gene_snATAC(ensemble, gene, grouping, outliers):
#def get_gene_snATAC(ensemble, gene, grouping, clustering, outliers, tsne_type='ATAC_ndim2_perp20'):
    """Return snATAC data points for a given gene.

    Data from ID-to-Name mapping and tSNE points are combined for plot generation.

    Arguments:
        ensemble (str): Name of ensemble.
        gene (str): Ensembl ID of gene.
        grouping (str): Variable for grouping cells. "cluster", "annotation", or "dataset".
        clustering (str): Different clustering algorithms and parameters. 'lv' = Louvain clustering.
        outliers (bool): Whether if outliers should be kept.
        tsne_type (str): Options for calculating tSNE. ndims = number of dimensions, perp = perplexity.

    Returns:
        DataFrame
    """

    # Prevent SQL injected since column names cannot be parameterized.
    if ";" in ensemble or ";" in grouping:
        return None

    # This query is just to fix gene id's missing the ensemble version number. 
    # Necessary because the table name must match exactly with whats on the MySQL database.
    # Ex. ENSMUSG00000026787 is fixed to ENSMUSG00000026787.3 -> gene_ENSMUSG00000026787_3 (table name in MySQL)
    result = db.get_engine(current_app, 'snATAC_data').execute("SELECT gene_id FROM genes WHERE gene_id LIKE %s", (gene+"%",)).fetchone()
    gene_table_name = 'gene_' + result['gene_id'].replace('.','_')
    
    query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, \
        %(ensemble)s.annotation_ATAC, %(ensemble)s.cluster_ATAC, \
        %(ensemble)s.tsne_x_ATAC, %(ensemble)s.tsne_y_ATAC, \
        %(gene_table_name)s.normalized_counts \
        FROM cells \
        INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
        LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble, 
                                                                                                'gene_table_name': gene_table_name,}
    # else:
    #     query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, cells.global_%(methylation_type)s, \
    #         %(ensemble)s.annotation_ATAC, %(ensemble)s.cluster_ATAC, \
    #         %(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, %(ensemble)s.tsne_z_%(tsne_type)s, \
    #         %(gene_table_name)s.normalized_counts \
    #         FROM cells \
    #         INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
    #         LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble, 
    #                                                                                                  'gene_table_name': gene_table_name,
    #                                                                                                  'tsne_type': tsne_type,
    #                                                                                                  'clustering': clustering,}
    try:
        df = pd.read_sql(query, db.get_engine(current_app, 'snATAC_data'))
    except exc.ProgrammingError as e:
        now = datetime.datetime.now()
        print("[{}] ERROR in app(get_gene_snATAC): {}".format(str(now), e))
        sys.stdout.flush()
        return None

    if df.empty: # If no data in column, return None 
        now = datetime.datetime.now()
        print("[{}] ERROR in app(get_gene_snATAC): No snATAC data for {}".format(str(now), ensemble))
        sys.stdout.flush()
        return None

    if grouping == 'annotation':
        df.fillna({'annotation_ATAC': 'None'}, inplace=True)
        df['annotation_cat'] = pd.Categorical(df['annotation_ATAC'], cluster_annotation_order)
        return df.sort_values(by='annotation_cat')
    elif grouping == 'cluster':
        return df.sort_values(by='cluster_ATAC')
    else:
        return df


@cache.memoize(timeout=1800)
def get_mult_gene_snATAC(ensemble, genes, grouping):
#def get_mult_gene_snATAC(ensemble, genes, grouping, clustering, tsne_type='ATAC_ndim2_perp20'):
    """Return averaged methylation data ponts for a set of genes.

    Data from ID-to-Name mapping and tSNE points are combined for plot generation.

    Arguments:
        ensemble (str): Name of ensemble.
        genes ([str]): List of gene IDs to query.
        grouping (str): Variable for grouping cells. "cluster", "annotation", or "dataset".
        clustering (str): Different clustering algorithms and parameters. 'lv' = Louvain clustering.
        tsne_type (str): Options for calculating tSNE. ndims = number of dimensions, perp = perplexity.

    Returns:
        DataFrame
    """

    # Prevent SQL injected since column names cannot be parameterized.
    if ";" in ensemble or ";" in grouping:
        return None

    genes = [gene+"%" for gene in genes]

    # This query is just to fix gene id's missing the ensemble version number. 
    # Necessary because the table name must match exactly with whats on the MySQL database.
    # Ex. ENSMUSG00000026787 is fixed to ENSMUSG00000026787.3
    first_query = "SELECT gene_id FROM genes WHERE gene_id LIKE %s" + " OR gene_id LIKE %s" * (len(genes)-1)
    result = db.get_engine(current_app, 'methylation_data').execute(first_query, (genes,)).fetchall()

    gene_table_names = ['gene_' + gene_id[0].replace('.','_') for gene_id in result]

    df_all = pd.DataFrame()
    
    first = True
    for gene_table_name in gene_table_names:
        query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, \
            %(ensemble)s.annotation_ATAC, %(ensemble)s.cluster_ATAC, \
            %(ensemble)s.tsne_x_ATAC, %(ensemble)s.tsne_y_ATAC, \
            %(gene_table_name)s.normalized_counts \
            FROM cells \
            INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
            LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble, 
                                                                                                    'gene_table_name': gene_table_name,}
        # else:
        #     query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, \
        #         %(ensemble)s.annotation_%(clustering)s, %(ensemble)s.cluster_%(clustering)s, \
        #         %(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, %(ensemble)s.tsne_z_%(tsne_type)s, \
        #         %(gene_table_name)s.normalized_counts \
        #         FROM cells \
        #         INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
        #         LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble, 
        #                                                                                                 'gene_table_name': gene_table_name,
        #                                                                                                 'tsne_type': tsne_type,
        #                                                                                                 'clustering': clustering,}
        try:
            df_all = df_all.append(pd.read_sql(query, db.get_engine(current_app, 'snATAC_data')))
        except exc.ProgrammingError as e:
            now = datetime.datetime.now()
            print("[{}] ERROR in app(get_mult_gene_snATAC): {}".format(str(now), e))
            sys.stdout.flush()
            return None
        
        if first:
            df_coords = df_all
        first = False

    if df_all.empty: # If no data in column, return None 
        now = datetime.datetime.now()
        print("[{}] ERROR in app(get_gene_snATAC): No snATAC data for {}".format(str(now), ensemble))
        sys.stdout.flush()
        return None

    df_avg_methylation = df_all.groupby(by='cell_id', as_index=False)['normalized_counts'].mean()
    df_coords.update(df_avg_methylation)

    if grouping == 'annotation':
        df_coords.fillna({'annotation_ATAC': 'None'}, inplace=True)
        df_coords['annotation_cat'] = pd.Categorical(df_coords['annotation_ATAC'], cluster_annotation_order)
        return df_coords.sort_values(by='annotation_cat')
    elif grouping == 'cluster':
        return df_coords.sort_values(by='cluster_ATAC')
    else:
        return df_coords


@cache.memoize(timeout=3600)
def get_ortholog_cluster_order():
    """Order cluster mm_hs_homologous_cluster.txt.

    Arguments:
        None

    Returns:
        list: tuples of (ensemble, cluster_number)
    """
    try:
        df = pd.read_csv(
            '{}/mm_hs_homologous_cluster.txt'.format(
                current_app.config['DATA_DIR'],
                engine='python'),
            sep='\t')
    except IOError:
        now = datetime.datetime.now()
        print("[{}] ERROR in app(get_ortholog_cluster_order): Could not load mm_hs_homologous_cluster.txt".format(str(now)))
        sys.stdout.flush()
        return []

    clusters = list()
    for _, row in df.iterrows():
        mmu_cluster = ('mmu', int(row['Mouse Cluster']))
        hsa_cluster = ('hsa', int(row['Human Cluster']))
        if mmu_cluster not in clusters:
            clusters.append(mmu_cluster)
        if hsa_cluster not in clusters:
            clusters.append(hsa_cluster)

    return clusters


@cache.memoize(timeout=1800)
def get_snATAC_scatter(ensemble, genes_query, grouping, ptile_start, ptile_end, tsne_outlier_bool):
#def get_snATAC_scatter(ensemble, tsne_type, genes_query, grouping, clustering, ptile_start, ptile_end, tsne_outlier_bool):
    """Generate scatter plot and gene body mCH scatter plot using tSNE coordinates from methylation(snmC-seq) data.

    Arguments:
        ensemble (str): Name of ensemble.
        tsne_type (str): Options for calculating tSNE. ndims = number of dimensions, perp = perplexity.
        genes_query (str):  Ensembl ID of gene(s) separated by spaces.
        grouping (str): Variable to group cells by. "cluster", "annotation", "dataset".
        clustering (str): Different clustering algorithms and parameters. 'lv' = Louvain clustering.
        ptile_start (float): Lower end of color percentile. [0, 1].
        ptile_end (float): Upper end of color percentile. [0, 1].
        tsne_outlier_bool (bool): Whether or not to change X and Y axes range to hide outliers. True = show outliers. 

    Returns:
        str: HTML generated by Plot.ly.
    """

    genes_query = genes_query.split()
    genes = []
    for gene in genes_query:
        genes.append(convert_gene_id_mmu_hsa(ensemble,gene))

    gene_name = ""
    x, y, text, mch = list(), list(), list(), list()

    if len(genes) == 1:
        points = get_gene_snATAC(ensemble, genes[0], grouping, True)
        gene_name = get_gene_by_id(ensemble, genes[0])['gene_name']
        title = 'Gene body snATAC counts: ' + gene_name
    else:
        points = get_mult_gene_snATAC(ensemble, genes, grouping)
        for gene in genes:
            gene_name += get_gene_by_id(ensemble, gene)['gene_name'] + '+'
        gene_name = gene_name[:-1]
        title = 'Avg. Gene body counts: ' + gene_name
    
    if points is None:
        raise FailToGraphException

    ### TSNE ### 
    if grouping != 'dataset':
        if grouping+'_ATAC' not in points.columns or points[grouping+'_ATAC'].nunique() <= 1:
            grouping = "cluster"
            print("**** Grouping by cluster")

    datasets = points['dataset'].unique().tolist()
    annotation_additional_y = 0.00 
    if grouping == 'dataset':
        unique_groups = datasets
        max_cluster = len(unique_groups)
    else:
        if grouping == 'cluster':
            annotation_additional_y = 0.025 # Necessary because legend items overlap with legend title (annotation) when there are many legend items
        max_cluster = points['cluster_ATAC'].max()
        unique_groups = points[grouping+'_ATAC'].unique().tolist()
    num_colors = len(unique_groups)
    
    colors = generate_cluster_colors(num_colors, grouping)
    symbols = ['circle', 'square', 'cross', 'triangle-up', 'triangle-down', 'octagon', 'star', 'diamond']
    
    traces_tsne = OrderedDict()

    legend_x = -.11
    grouping_clustering = grouping
    if grouping != 'dataset':
        legend_x = -.10
        grouping_clustering = grouping+'_ATAC'

    unique_groups = points[grouping_clustering].unique().tolist()

    if tsne_outlier_bool:
        top_x = points['tsne_x_ATAC'].quantile(0.999)
        bottom_x = points['tsne_x_ATAC'].quantile(0.001) 
        top_y = points['tsne_y_ATAC'].quantile(0.999)
        bottom_y = points['tsne_y_ATAC'].quantile(0.001) 
        
    else:
        top_x = points['tsne_x_ATAC'].max()
        bottom_x = points['tsne_x_ATAC'].min()
        top_y = points['tsne_y_ATAC'].max()
        bottom_y = points['tsne_y_ATAC'].min()

    range_x = top_x - bottom_x 
    top_x = top_x + range_x * 0.1
    bottom_x = bottom_x - range_x*0.1
    range_y = top_y - bottom_y
    top_y = top_y + range_y * 0.1
    bottom_y = bottom_y - range_y*0.1

    if len(points) > 3000: 
        marker_size = 2
    else:
        marker_size = 4

    ## 2D tSNE coordinates ##
    for i, group in enumerate(unique_groups):
        points_group = points[points[grouping_clustering]==group]
        if grouping_clustering.startswith('cluster'):
            group_str = 'cluster_' + str(group)
        elif grouping_clustering== "dataset":
            group = group.strip('CEMBA_')
            group_str = group
        else:
            group_str = group

        color_num = i
        
        trace2d = traces_tsne.setdefault(color_num, Scatter(
            x=list(),
            y=list(),
            text=list(),
            mode='markers',
            visible=True,
            name=group_str,
            legendgroup=group,
            marker={
                   'color': colors[color_num],
                   'size': marker_size,
                   #'opacity': 0.8,
                   #'symbol': symbols[datasets.index(dataset)],
            },
            hoverinfo='text'))
        trace2d['x'] = points_group['tsne_x_ATAC'].values.tolist()
        trace2d['y'] = points_group['tsne_y_ATAC'].values.tolist()
        trace2d['text'] = [build_hover_text({'Cell Name': point[1],
                                             'Dataset': point[2],
                                             'Annotation': point[3],
                                             'Cluster': point[4]})
                           for point in points_group.itertuples(index=False)]

    ### snATAC reads SCATTER ### 
    x = points['tsne_x_ATAC'].tolist()
    y = points['tsne_y_ATAC'].tolist()
    ATAC_counts = points['normalized_counts']
    text_ATAC = [build_hover_text({'Normalized Counts': round(point[-1], 6),
                                   'Cell Name': point[1],
                                   'Annotation': point[3],
                                   'Cluster': point[4]})
                 for point in points.itertuples(index=False)]


    ATAC_dataframe = pd.DataFrame(ATAC_counts)
    start = ATAC_dataframe.dropna().quantile(ptile_start)[0].tolist()
    end = ATAC_dataframe.dropna().quantile(ptile_end).values[0].tolist()
    ATAC_colors = [set_color_by_percentile(x, start, end) for x in ATAC_counts]

    colorbar_tickval = list(arange(start, end, (end - start) / 4))
    colorbar_tickval[0] = start
    colorbar_tickval.append(end)
    colorbar_ticktext = [
        str(round(x, 3)) for x in arange(start, end, (end - start) / 4)
    ]
    colorbar_ticktext[0] = '<' + str(round(start, 3))
    colorbar_ticktext.append('>' + str(round(end, 3)))

    trace_ATAC = Scatter(
        mode='markers',
        x=x,
        y=y,
        text=text_ATAC,
        marker={
            'color': ATAC_colors,
            'colorscale': 'Viridis',
            'size': marker_size,
            'colorbar': {
                'x': 1.05,
                'len': 0.5,
                'thickness': 10,
                'title': 'Normalized Counts',
                'titleside': 'right',
                'tickmode': 'array',
                'tickvals': colorbar_tickval,
                'ticktext': colorbar_ticktext,
                'tickfont': {'size': 10}
            }
        },
        showlegend=False,
        yaxis='y',
        xaxis='x2',
        hoverinfo='text')

    layout = Layout(
        autosize=True,
        height=450,
        title=title,
        titlefont={'color': 'rgba(1,2,2,1)',
                   'size': 16},
        legend={'x':legend_x,
                'y':0.95,
                'tracegroupgap': 0.5},
        margin={'l': 49,
                'r': 0,
                'b': 30,
                't': 50,
                'pad': 0},
        xaxis={
            'domain': [0, 0.49],
            'type': 'linear',
            'ticks': '',
            'dtick': 10,
            'tickwidth': 0,
            'showticklabels': False,
            'showline': False,
            'showgrid': True,
            'zeroline': False,
            'linecolor': 'black',
            'linewidth': 0.5,
            'mirror': False,
            'scaleanchor': 'x2',
            'range':[bottom_x,top_x]
        },
        xaxis2={
            'domain': [0.51, 1],
            'type': 'linear',
            'ticks': '',
            'dtick': 10,
            'tickwidth': 0,
            'showticklabels': False,
            'showline': False,
            'showgrid': True,
            'zeroline': False,
            'linecolor': 'black',
            'linewidth': 0.5,
            'mirror': False,
            'scaleanchor': 'y',
            'range':[bottom_x,top_x]
        },
        yaxis={
            'domain': [0,1],
            'type': 'linear',
            'ticks': '',
            'dtick': 10,
            'tickwidth': 0,
            'showticklabels': False,
            'showline': False,
            'showgrid': True,
            'side': 'right',
            'zeroline': False,
            'linecolor': 'black',
            'linewidth': 0.5,
            'mirror': False,
            'range':[bottom_y,top_y]
        },
        hovermode='closest',
        hoverdistance='10',)

    
    fig = tools.make_subplots(
            rows=1,
            cols=2,
            shared_xaxes=False,
            shared_yaxes=True,
            print_grid=False,
            subplot_titles=("tSNE", "Normalized Counts"),
            )

    for trace in traces_tsne.items():
        fig.append_trace(trace[1], 1,1)
    fig.append_trace(trace_ATAC, 1,2)

    fig['layout'].update(layout)
    fig['layout']['annotations'].extend([Annotation(text="Cluster Labels",
                                                    x=legend_x+0.01,
                                                    y=1.02 + annotation_additional_y,
                                                    xanchor="left",
                                                    yanchor="top",
                                                    showarrow=False,
                                                    xref="paper",
                                                    yref="paper",
                                                    font={'size': 12,
                                                          'color': 'gray',})])



    ## 3D tSNE coordinates ##
    # else: 
    #     for i, group in enumerate(unique_groups):

    #         points_group = points[points[grouping_clustering]==group]
    #         if grouping_clustering.startswith('cluster'):
    #             group_str = 'cluster_' + str(group)
    #         elif grouping_clustering== "dataset":
    #             group = group.strip('CEMBA_')
    #             group_str = group
    #         else:
    #             group_str = group

    #         color_num = i
    #         
    #         trace3d = traces_tsne.setdefault(color_num, Scatter3d(
    #             x=list(),
    #             y=list(),
    #             z=list(),
    #             text=list(),
    #             mode='markers',
    #             visible=True,
    #             name=group_str,
    #             legendgroup=group,
    #             scene='scene1',
    #             marker={
    #                    'color': colors[color_num],
    #                    'size': marker_size,
    #                    'opacity': 0.8,
    #                    #'symbol': symbols[datasets.index(dataset)],
    #             },
    #             hoverinfo='text'))
    #         trace3d['x'] = points_group['tsne_x_'+tsne_type].values.tolist()
    #         trace3d['y'] = points_group['tsne_y_'+tsne_type].values.tolist()
    #         trace3d['z'] = points_group['tsne_z_'+tsne_type].values.tolist()
    #         trace3d['text'] = [build_hover_text({'Cell Name': point[1],
    #                                              'Dataset': point[2],
    #                                              'Annotation': point[3],
    #                                              'Cluster': point[4]})
    #                            for point in points_group.itertuples(index=False)]

    #     ### sNATAC COUNTS SCATTER ### 
    #     x = points['tsne_x_' + tsne_type].tolist()
    #     y = points['tsne_y_' + tsne_type].tolist()
    #     z = points['tsne_z_' + tsne_type].tolist()
    #     ATAC_counts = points['normalized_counts']
    #     text_ATAC = [build_hover_text({'Normalized Counts': round(point[-1], 6),
    #                                    'Cell Name': point[1],
    #                                    'Annotation': point[3],
    #                                    'Cluster': point[4]})
    #                  for point in points.itertuples(index=False)]


    #     ATAC_dataframe = pd.DataFrame(ATAC_counts)
    #     start = ATAC_dataframe.dropna().quantile(ptile_start)[0].tolist()
    #     end = ATAC_dataframe.dropna().quantile(ptile_end).values[0].tolist()
    #     ATAC_colors = [set_color_by_percentile(x, start, end) for x in ATAC_counts]

    #     colorbar_tickval = list(arange(start, end, (end - start) / 4))
    #     colorbar_tickval[0] = start
    #     colorbar_tickval.append(end)
    #     colorbar_ticktext = [
    #         str(round(x, 3)) for x in arange(start, end, (end - start) / 4)
    #     ]
    #     colorbar_ticktext[0] = '<' + str(round(start, 3))
    #     colorbar_ticktext.append('>' + str(round(end, 3)))

    #     trace_ATAC = Scatter3d(
    #         mode='markers',
    #         x=x,
    #         y=y,
    #         z=z,
    #         text=text_ATAC,
    #         scene='scene2',
    #         marker={
    #             'color': ATAC_colors,
    #             'colorscale': 'Viridis',
    #             'size': marker_size,
    #             'colorbar': {
    #                 'x': 1.05,
    #                 'len': 0.5,
    #                 'thickness': 10,
    #                 'title': 'Normalized Counts',
    #                 'titleside': 'right',
    #                 'tickmode': 'array',
    #                 'tickvals': colorbar_tickval,
    #                 'ticktext': colorbar_ticktext,
    #                 'tickfont': {'size': 10}
    #             }
    #         },
    #         showlegend=False,
    #         hoverinfo='text')

    #     layout = Layout(
    #         autosize=True,
    #         height=450,
    #         title=title,
    #         titlefont={'color': 'rgba(1,2,2,1)',
    #                    'size': 16},
    #         legend={'x':-.1, 'y':1},
    #         margin={'l': 49,
    #                 'r': 0,
    #                 'b': 30,
    #                 't': 50,
    #                 'pad': 10
    #                 },

    #         hovermode='closest',
    #     )

    #     scene={
    #         'camera':{
    #             'eye': dict(x=1.2, y=1.5, z=0.7),
    #             'center': dict(x=0.25, z=-0.1)
    #                  },
    #         'aspectmode':'data',
    #         'xaxis':{
    #             'title': 'tSNE 1',
    #             'titlefont': {
    #                 'color': 'rgba(1,2,2,1)',
    #                 'size': 12
    #             },
    #             'type': 'linear',
    #             'ticks': '',
    #             'showticklabels': False,
    #             'tickwidth': 0,
    #             'showline': True,
    #             'showgrid': False,
    #             'zeroline': False,
    #             'linecolor': 'black',
    #             'linewidth': 0.5,
    #             'mirror': True
    #         },
    #         'yaxis':{
    #             'title': 'tSNE 2',
    #             'titlefont': {
    #                 'color': 'rgba(1,2,2,1)',
    #                 'size': 12
    #             },
    #             'type': 'linear',
    #             'ticks': '',
    #             'showticklabels': False,
    #             'tickwidth': 0,
    #             'showline': True,
    #             'showgrid': False,
    #             'zeroline': False,
    #             'linecolor': 'black',
    #             'linewidth': 0.5,
    #             'mirror': True
    #         },
    #         'zaxis':{
    #             'title': 'tSNE 3',
    #             'titlefont': {
    #                 'color': 'rgba(1,2,2,1)',
    #                 'size': 12
    #             },
    #             'type': 'linear',
    #             'ticks': '',
    #             'showticklabels': False,
    #             'tickwidth': 0,
    #             'showline': True,
    #             'showgrid': False,
    #             'zeroline': False,
    #             'linecolor': 'black',
    #             'linewidth': 0.5,
    #             'mirror': True
    #         },
    #     }

    #     fig = tools.make_subplots(rows=1,
    #                               cols=2,
    #                               shared_xaxes=True,
    #                               #shared_yaxes=True,
    #                               print_grid=False,
    #                               subplot_titles=("tSNE", "Normalized_Counts"),
    #                               specs=[[{'is_3d':True}, {'is_3d':True}]])

    #     for trace in traces_tsne.items():
    #         fig.append_trace(trace[1], 1,1)
    #     fig.append_trace(trace_ATAC, 1,2)

    #     fig['layout'].update(layout)
    #     fig['layout']['scene1'].update(scene)
    #     fig['layout']['scene2'].update(scene)
    # 
    #     fig['layout']['annotations'].extend([Annotation(text="Cluster Labels",
    #                                                     x=-.09,
    #                                                     y=1.03 + annotation_additional_y,
    #                                                     xanchor="left",
    #                                                     yanchor="top",
    #                                                     showarrow=False,
    #                                                     xref="paper",
    #                                                     yref="paper",
    #                                                     font={'size': 12,
    #                                                           'color': 'gray',})])

    return plotly.offline.plot(
        figure_or_data=fig,
        output_type='div',
        show_link=False,
        include_plotlyjs=False)


@cache.memoize(timeout=1800)
def get_methylation_scatter(ensemble, tsne_type, methylation_type, genes_query, level, grouping, clustering, ptile_start, ptile_end, tsne_outlier_bool):
    """Generate scatter plot and gene body reads scatter plot using tSNE coordinates from snATAC-seq data.

    Arguments:
        ensemble (str): Name of ensemble.
        tsne_type (str): Options for calculating tSNE. ndims = number of dimensions, perp = perplexity.
        methylation_type (str): Type of methylation to visualize. "mCH", "mCG", or "mCA".
        genes_query (str):  Ensembl ID of gene(s) separated by spaces.
        level (str): "original" or "normalized" methylation values.
        grouping (str): Variable to group cells by. "cluster", "annotation", "dataset".
        clustering (str): Different clustering algorithms and parameters. 'lv' = Louvain clustering.
        ptile_start (float): Lower end of color percentile. [0, 1].
        ptile_end (float): Upper end of color percentile. [0, 1].
        tsne_outlier_bool (bool): Whether or not to change X and Y axes range to hide outliers. True = do show outliers. 

    Returns:
        str: HTML generated by Plot.ly.
    """


    genes_query = genes_query.split()
    genes = []
    context = methylation_type[1:]
    for gene in genes_query:
        genes.append(convert_gene_id_mmu_hsa(ensemble,gene))

    gene_name = ""
    x, y, text, mch = list(), list(), list(), list()

    if len(genes) == 1:
        points = get_gene_methylation(ensemble, methylation_type, genes[0], grouping, clustering, level, True, tsne_type)
        gene_name = get_gene_by_id(ensemble, genes[0])['gene_name']
        title = 'Gene body ' + methylation_type + ': ' + gene_name
    else:
        points = get_mult_gene_methylation(ensemble, methylation_type, genes, grouping, clustering, level, tsne_type)
        for gene in genes:
            gene_name += get_gene_by_id(ensemble, gene)['gene_name'] + '+'
        gene_name = gene_name[:-1]
        title = 'Avg. Gene body ' + methylation_type + ': ' + gene_name

    if points is None:
        raise FailToGraphException

    ### TSNE ### 
    if grouping != 'dataset':
        if grouping+'_'+clustering not in points.columns or points[grouping+'_'+clustering].nunique() <= 1:
            grouping = "cluster"
            clustering = "mCH_lv_npc50_k30"
            print("**** Using cluster_mCH_lv_npc50_k30")

    datasets = points['dataset'].unique().tolist()
    annotation_additional_y = 0.00 
    if grouping == 'dataset':
        unique_groups = datasets
        max_cluster = len(unique_groups)
    else:
        if grouping == 'cluster':
            annotation_additional_y = 0.025 # Necessary because legend items overlap with legend title (annotation) when there are many legend items
        max_cluster = points['cluster_'+clustering].max()
        unique_groups = points[grouping+'_'+clustering].unique().tolist()
    num_colors = len(unique_groups)
    
    colors = generate_cluster_colors(num_colors, grouping)
    symbols = ['circle', 'square', 'cross', 'triangle-up', 'triangle-down', 'octagon', 'star', 'diamond']
    
    traces_tsne = OrderedDict()

    legend_x = -.11
    grouping_clustering = grouping
    if grouping != 'dataset':
        legend_x = -.10
        grouping_clustering = grouping+'_'+clustering

    unique_groups = points[grouping_clustering].unique().tolist()

    if tsne_outlier_bool:
        top_x = points['tsne_x_'+tsne_type].quantile(0.999)
        bottom_x = points['tsne_x_'+tsne_type].quantile(0.001) 
        top_y = points['tsne_y_'+tsne_type].quantile(0.999)
        bottom_y = points['tsne_y_'+tsne_type].quantile(0.001) 
        
    else:
        top_x = points['tsne_x_'+tsne_type].max()
        bottom_x = points['tsne_x_'+tsne_type].min()
        top_y = points['tsne_y_'+tsne_type].max()
        bottom_y = points['tsne_y_'+tsne_type].min()

    range_x = top_x - bottom_x 
    top_x = top_x + range_x * 0.1
    bottom_x = bottom_x - range_x*0.1
    range_y = top_y - bottom_y
    top_y = top_y + range_y * 0.1
    bottom_y = bottom_y - range_y*0.1

    if len(points) > 3000: 
        marker_size = 2
    else:
        marker_size = 4

    ## 2D tSNE coordinates ##
    if 'ndim2' in tsne_type:
        for i, group in enumerate(unique_groups):
            points_group = points[points[grouping_clustering]==group]
            if grouping_clustering.startswith('cluster'):
                group_str = 'cluster_' + str(group)
            elif grouping_clustering== "dataset":
                group = group.strip('CEMBA_')
                group_str = group
            else:
                group_str = group

            color_num = i
            
            trace2d = traces_tsne.setdefault(color_num, Scatter(
                x=list(),
                y=list(),
                text=list(),
                mode='markers',
                visible=True,
                name=group_str,
                legendgroup=group,
                marker={
                       'color': colors[color_num],
                       'size': marker_size,
                       #'opacity': 0.8,
                       #'symbol': symbols[datasets.index(dataset)],
                },
                hoverinfo='text'))
            trace2d['x'] = points_group['tsne_x_'+tsne_type].values.tolist()
            trace2d['y'] = points_group['tsne_y_'+tsne_type].values.tolist()
            trace2d['text'] = [build_hover_text({'Cell Name': point[1],
                                                 'Dataset': point[2],
                                                 'Annotation': point[4],
                                                 'Cluster': point[5]})
                               for point in points_group.itertuples(index=False)]

        ### METHYLATION SCATTER ### 
        x = points['tsne_x_' + tsne_type].tolist()
        y = points['tsne_y_' + tsne_type].tolist()
        mch = points[methylation_type + '/' + context + '_' + level]
        text_methylation = [build_hover_text({level.title()+' '+methylation_type: round(point[-2], 6),
                                  'Cell Name': point[1],
                                  'Annotation': point[4],
                                  'Cluster': point[5]})
                for point in points.itertuples(index=False)]


        mch_dataframe = pd.DataFrame(mch)
        start = mch_dataframe.dropna().quantile(ptile_start)[0].tolist()
        end = mch_dataframe.dropna().quantile(ptile_end).values[0].tolist()
        mch_colors = [set_color_by_percentile(x, start, end) for x in mch]

        colorbar_tickval = list(arange(start, end, (end - start) / 4))
        colorbar_tickval[0] = start
        colorbar_tickval.append(end)
        colorbar_ticktext = [
            str(round(x, 3)) for x in arange(start, end, (end - start) / 4)
        ]
        colorbar_ticktext[0] = '<' + str(round(start, 3))
        colorbar_ticktext.append('>' + str(round(end, 3)))

        trace_methylation = Scatter(
            mode='markers',
            x=x,
            y=y,
            text=text_methylation,
            marker={
                'color': mch_colors,
                'colorscale': 'Viridis',
                'size': marker_size,
                'colorbar': {
                    'x': 1.05,
                    'len': 0.5,
                    'thickness': 10,
                    'title': level.capitalize() + ' ' + methylation_type,
                    'titleside': 'right',
                    'tickmode': 'array',
                    'tickvals': colorbar_tickval,
                    'ticktext': colorbar_ticktext,
                    'tickfont': {'size': 10}
                }
            },
            showlegend=False,
            yaxis='y',
            xaxis='x2',
            hoverinfo='text')

        layout = Layout(
            autosize=True,
            height=450,
            title=title,
            titlefont={'color': 'rgba(1,2,2,1)',
                       'size': 16},
            legend={'x':legend_x,
                    'y':0.95,
                    'tracegroupgap': 0.5},
            margin={'l': 49,
                    'r': 0,
                    'b': 30,
                    't': 50,
                    'pad': 0},
            xaxis={
                'domain': [0, 0.49],
                'type': 'linear',
                'ticks': '',
                'dtick': 10,
                'tickwidth': 0,
                'showticklabels': False,
                'showline': False,
                'showgrid': True,
                'zeroline': False,
                'linecolor': 'black',
                'linewidth': 0.5,
                'mirror': False,
                'scaleanchor': 'x2',
                'range':[bottom_x,top_x]
            },
            xaxis2={
                'domain': [0.51, 1],
                'type': 'linear',
                'ticks': '',
                'dtick': 10,
                'tickwidth': 0,
                'showticklabels': False,
                'showline': False,
                'showgrid': True,
                'zeroline': False,
                'linecolor': 'black',
                'linewidth': 0.5,
                'mirror': False,
                'scaleanchor': 'y',
                'range':[bottom_x,top_x]
            },
            yaxis={
                'domain': [0,1],
                'type': 'linear',
                'ticks': '',
                'dtick': 10,
                'tickwidth': 0,
                'showticklabels': False,
                'showline': False,
                'showgrid': True,
                'side': 'right',
                'zeroline': False,
                'linecolor': 'black',
                'linewidth': 0.5,
                'mirror': False,
                'range':[bottom_y,top_y]
            },
            hovermode='closest',
            hoverdistance='10',)

        
        fig = tools.make_subplots(
                rows=1,
                cols=2,
                shared_xaxes=False,
                shared_yaxes=True,
                print_grid=False,
                subplot_titles=("tSNE", level.title()+" Methylation ("+methylation_type+")"),
                )

        for trace in traces_tsne.items():
            fig.append_trace(trace[1], 1,1)
        fig.append_trace(trace_methylation, 1,2)

        fig['layout'].update(layout)
        fig['layout']['annotations'].extend([Annotation(text="Cluster Labels",
                                                        x=legend_x+0.01,
                                                        y=1.02 + annotation_additional_y,
                                                        xanchor="left",
                                                        yanchor="top",
                                                        showarrow=False,
                                                        xref="paper",
                                                        yref="paper",
                                                        font={'size': 12,
                                                              'color': 'gray',})])



    ## 3D tSNE coordinates ##
    else: 
        for i, group in enumerate(unique_groups):

            points_group = points[points[grouping_clustering]==group]
            if grouping_clustering.startswith('cluster'):
                group_str = 'cluster_' + str(group)
            elif grouping_clustering== "dataset":
                group = group.strip('CEMBA_')
                group_str = group
            else:
                group_str = group

            color_num = i
            
            trace3d = traces_tsne.setdefault(color_num, Scatter3d(
                x=list(),
                y=list(),
                z=list(),
                text=list(),
                mode='markers',
                visible=True,
                name=group_str,
                legendgroup=group,
                scene='scene1',
                marker={
                       'color': colors[color_num],
                       'size': marker_size,
                       'opacity': 0.8,
                       #'symbol': symbols[datasets.index(dataset)],
                },
                hoverinfo='text'))
            trace3d['x'] = points_group['tsne_x_'+tsne_type].values.tolist()
            trace3d['y'] = points_group['tsne_y_'+tsne_type].values.tolist()
            trace3d['z'] = points_group['tsne_z_'+tsne_type].values.tolist()
            trace3d['text'] = [build_hover_text({'Cell Name': point[1],
                                                 'Dataset': point[2],
                                                 'Annotation': point[4],
                                                 'Cluster': point[5]})
                               for point in points_group.itertuples(index=False)]

        ### METHYLATION SCATTER ### 
        x = points['tsne_x_' + tsne_type].tolist()
        y = points['tsne_y_' + tsne_type].tolist()
        z = points['tsne_z_' + tsne_type].tolist()
        mch = points[methylation_type + '/' + context + '_' + level]
        text_methylation = [build_hover_text({methylation_type: round(point[-2], 6),
                                  'Cell Name': point[1],
                                  'Annotation': point[4],
                                  'Cluster': point[5]})
                for point in points.itertuples(index=False)]


        mch_dataframe = pd.DataFrame(mch)
        start = mch_dataframe.dropna().quantile(ptile_start)[0].tolist()
        end = mch_dataframe.dropna().quantile(ptile_end).values[0].tolist()
        mch_colors = [set_color_by_percentile(x, start, end) for x in mch]

        colorbar_tickval = list(arange(start, end, (end - start) / 4))
        colorbar_tickval[0] = start
        colorbar_tickval.append(end)
        colorbar_ticktext = [
            str(round(x, 3)) for x in arange(start, end, (end - start) / 4)
        ]
        colorbar_ticktext[0] = '<' + str(round(start, 3))
        colorbar_ticktext.append('>' + str(round(end, 3)))

        trace_methylation = Scatter3d(
            mode='markers',
            x=x,
            y=y,
            z=z,
            text=text_methylation,
            scene='scene2',
            marker={
                'color': mch_colors,
                'colorscale': 'Viridis',
                'size': marker_size,
                'colorbar': {
                    'x': 1.05,
                    'len': 0.5,
                    'thickness': 10,
                    'title': level.capitalize() + ' ' + methylation_type,
                    'titleside': 'right',
                    'tickmode': 'array',
                    'tickvals': colorbar_tickval,
                    'ticktext': colorbar_ticktext,
                    'tickfont': {'size': 10}
                }
            },
            showlegend=False,
            hoverinfo='text')

        layout = Layout(
            autosize=True,
            height=450,
            title=title,
            titlefont={'color': 'rgba(1,2,2,1)',
                       'size': 16},
            legend={'x':-.1, 'y':1},
            margin={'l': 49,
                    'r': 0,
                    'b': 30,
                    't': 50,
                    'pad': 10
                    },

            hovermode='closest',
        )

        scene={
            'camera':{
                'eye': dict(x=1.2, y=1.5, z=0.7),
                'center': dict(x=0.25, z=-0.1)
                     },
            'aspectmode':'data',
            'xaxis':{
                'title': 'tSNE 1',
                'titlefont': {
                    'color': 'rgba(1,2,2,1)',
                    'size': 12
                },
                'type': 'linear',
                'ticks': '',
                'showticklabels': False,
                'tickwidth': 0,
                'showline': True,
                'showgrid': False,
                'zeroline': False,
                'linecolor': 'black',
                'linewidth': 0.5,
                'mirror': True
            },
            'yaxis':{
                'title': 'tSNE 2',
                'titlefont': {
                    'color': 'rgba(1,2,2,1)',
                    'size': 12
                },
                'type': 'linear',
                'ticks': '',
                'showticklabels': False,
                'tickwidth': 0,
                'showline': True,
                'showgrid': False,
                'zeroline': False,
                'linecolor': 'black',
                'linewidth': 0.5,
                'mirror': True
            },
            'zaxis':{
                'title': 'tSNE 3',
                'titlefont': {
                    'color': 'rgba(1,2,2,1)',
                    'size': 12
                },
                'type': 'linear',
                'ticks': '',
                'showticklabels': False,
                'tickwidth': 0,
                'showline': True,
                'showgrid': False,
                'zeroline': False,
                'linecolor': 'black',
                'linewidth': 0.5,
                'mirror': True
            },
        }

        fig = tools.make_subplots(rows=1,
                                  cols=2,
                                  shared_xaxes=True,
                                  #shared_yaxes=True,
                                  print_grid=False,
                                  subplot_titles=("tSNE", "Methylation"),
                                  specs=[[{'is_3d':True}, {'is_3d':True}]])

        for trace in traces_tsne.items():
            fig.append_trace(trace[1], 1,1)
        fig.append_trace(trace_methylation, 1,2)

        fig['layout'].update(layout)
        fig['layout']['scene1'].update(scene)
        fig['layout']['scene2'].update(scene)
    
        fig['layout']['annotations'].extend([Annotation(text="Cluster Labels",
                                                        x=-.09,
                                                        y=1.03 + annotation_additional_y,
                                                        xanchor="left",
                                                        yanchor="top",
                                                        showarrow=False,
                                                        xref="paper",
                                                        yref="paper",
                                                        font={'size': 12,
                                                              'color': 'gray',})])

    return plotly.offline.plot(
        figure_or_data=fig,
        output_type='div',
        show_link=False,
        include_plotlyjs=False)


@cache.memoize(timeout=3600)
def get_mch_heatmap(ensemble, methylation_type, grouping, clustering, level, ptile_start, ptile_end, normalize_row, query):
    """Generate mCH heatmap comparing multiple genes.

    Arguments:
        ensemble (str): Name of ensemble.
        methylation_type (str): Type of methylation to visualize.        "mch" or "mcg"
        level (str): "original" or "normalized" methylation values.
        outliers (bool): Whether if outliers should be displayed.
        ptile_start (float): Lower end of color percentile. [0, 1].
        ptile_end (float): Upper end of color percentile. [0, 1].
        normalize_row (bool): Whether to normalize by each row (gene). 
        query ([str]): Ensembl IDs of genes to display.

    Returns:
        str: HTML generated by Plot.ly
    """

    tsne_type = 'mCH_ndim2_perp20'

    if normalize_row:
        normal_or_original = 'Normalized'
    else:
        normal_or_original = 'Original'

    title = normal_or_original + " gene body " + methylation_type + " by cluster: "
    genes = [convert_gene_id_mmu_hsa(ensemble,gene) for gene in query.split()]

    gene_info_df = pd.DataFrame()
    for gene_id in genes:
        gene_name = get_gene_by_id(ensemble, gene_id)['gene_name']
        title += gene_name + "+"
        gene_info_df[gene_name] = median_cluster_mch(get_gene_methylation(ensemble, methylation_type, gene_id, grouping, clustering, level, True), grouping, clustering)

    title = title[:-1] # Gets rid of last '+'

    gene_info_df.reset_index(inplace=True)
    if grouping == 'annotation':
        gene_info_df['annotation_cat'] = pd.Categorical(gene_info_df['annotation_'+clustering], cluster_annotation_order)
        gene_info_df.sort_values(by='annotation_cat', inplace=True)
        gene_info_df.drop('annotation_cat', axis=1, inplace=True)
    elif grouping == 'cluster':
        gene_info_df.sort_values(by='cluster_'+clustering, inplace=True)
    gene_info_df.set_index(grouping+'_'+clustering, inplace=True)

    normal_or_original = 'Original'
    if normalize_row:
        for gene in gene_info_df:
            # z-score
            # gene_info_df[gene] = (gene_info_df[gene] - gene_info_df[gene].mean()) / gene_info_df[gene].std()
            # min-max
            gene_info_df[gene] = (gene_info_df[gene] - gene_info_df[gene].min()) / (gene_info_df[gene].max() - gene_info_df[gene].min())
        normal_or_original = 'Normalized'

    gene_info_dict = gene_info_df.to_dict(into=OrderedDict)    

    x, y, text, hover, mch = list(), list(), list(), list(), list()
    i = 0
    name_prepend = ""
    if grouping == 'cluster':
        name_prepend = 'cluster_'
    for key in list(gene_info_dict.keys()):
        j = 0
        y.append(key)
        mch.append(list(gene_info_dict[key].values()))
        for cluster in list(gene_info_dict[key].keys()):
            x.append(name_prepend+str(cluster))
            text.append(build_hover_text({
                'Gene': key,
                'Cluster': x[j],
                methylation_type: mch[i][j]
            }))
            j += 1
        hover.append(text)
        text = []
        i += 1

    flat_mch = list(chain.from_iterable(mch))
    mch_dataframe = pd.DataFrame(flat_mch).dropna()
    start = mch_dataframe.quantile(0.05)[0].tolist()
    end = mch_dataframe.quantile(0.95).values[0].tolist()

    colorbar_tickval = list(arange(start, end, (end - start) / 4))
    colorbar_tickval[0] = start
    colorbar_tickval.append(end)
    colorbar_ticktext = [
        str(round(x, 3)) for x in arange(start, end, (end - start) / 4)
    ]
    if normalize_row == True:
        colorbar_ticktext[0] = str(round(start, 3))
    else:
        colorbar_ticktext[0] = '<' + str(round(start, 3))
    colorbar_ticktext.append('>' + str(round(end, 3)))

    # Due to a weird bug(?) in plotly, the number of elements in tickvals and ticktext 
    # must be greater than or equal to number of genes in query. Else, javascript throws 
    # Uncaught Typeerrors when trying to hover over genes. (Tomo 12/11/17)
    while len(colorbar_tickval) < len(genes):
        colorbar_tickval.insert(0,start)
        if normalize_row == True:
            colorbar_ticktext.insert(0, str(round(start, 3)))
        else:
            colorbar_ticktext.insert(0, '<' + str(round(start, 3)))

    trace = Heatmap(
        x=x,
        y=y,
        z=mch,
        text=hover,
        colorscale='Viridis',
        colorbar={
            'x': 1.0,
            'len': 1,
            'title': level.capitalize() + ' ' + methylation_type,
            'titleside': 'right',
            'tickmode': 'array',
            'tickvals': colorbar_tickval,
            'ticktext': colorbar_ticktext,
            'thickness': 10,
            'tickfont': {'size': 10}
            },
        hoverinfo='text'
        )

    layout = Layout(
        title=title,
        height=450,
        titlefont={'color': 'rgba(1,2,2,1)',
                   'size': 16},
        autosize=True,
        xaxis={
            'side': 'bottom',
            'tickangle': -45,
            'tickfont': {'size': 12}
               },
        yaxis={
            'title': 'Genes',
            'tickangle': 15,
            'tickfont': {'size': 12}
            },
        hovermode='closest'
        )

    # Available colorscales:
    # https://community.plot.ly/t/what-colorscales-are-available-in-plotly-and-which-are-the-default/2079
    updatemenus = list([
        dict(
            buttons=list([
                dict(
                    args=['colorscale', 'Viridis'],
                    label='Viridis',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Bluered'],
                    label='Bluered',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Blackbody'],
                    label='Blackbody',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Electric'],
                    label='Electric',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Earth'],
                    label='Earth',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Jet'],
                    label='Jet',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Rainbow'],
                    label='Rainbow',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Picnic'],
                    label='Picnic',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Portland'],
                    label='Portland',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'YlGnBu'],
                    label='YlGnBu',
                    method='restyle'
                )
            ]),
            direction='down',
            showactive=True,
            x=0,
            xanchor='left',
            y=1.17,
            yanchor='top'
        )
    ])

    layout['updatemenus'] = updatemenus

    return plotly.offline.plot(
        {
            'data': [trace],
            'layout': layout
        },
        output_type='div',
        show_link=False,
        include_plotlyjs=False)


def get_mch_heatmap_two_ensemble(ensemble, methylation_type, level, ptile_start, ptile_end, normalize_row, query):
    """Generate gene body mCH heatmap for two ensemble.

    Traces are grouped by cluster and ordered by mm_hs_homologous_cluster.txt.

    Arguments:
        ensemble (str): Species being viewed.
        methylation_type (str): Type of methylation to visualize. "mch" or "mcg"
        level (str): "original" or "normalized" methylation values.
        ptile_start (float): Lower end of color percentile. [0, 1].
        ptile_end (float): Upper end of color percentile. [0, 1].
        normalize_row (bool): Whether to normalize by each row (gene). 
        query ([str]):  List of ensembl IDs of genes.

    Returns:
        str: HTML generated by Plot.ly.
    """
    
    """ NOTES
    
        Must be on same colorscale, which means original heatmap must be replotted
            Color bar should be to right of right most map. Subplots?
                https://plot.ly/python/subplots/
        Ideally two separated 
        Left/Right vs. Up/Down
        if ortholog doesn't exist, fill with NaN
        Titles?
        Y-axis?

    """

    gene_mch_hsa_df = pd.DataFrame()
    gene_mch_mmu_df = pd.DataFrame()
    gene_mch_combined_df = pd.DataFrame()

    if ensemble == 'mouse_published' or ensemble == 'mmu':
        gene_id_list_mmu = [ convert_gene_id_mmu_hsa(ensemble, gene_id) for gene_id in query.split() ]
        gene_id_list_hsa = [ find_orthologs(mmu_gene_id = gene_id)['hsa_gene_id'] for gene_id in gene_id_list_mmu ]
        for gene_id in gene_id_list_mmu:
            gene_label_mmu = get_gene_by_id('mouse_published', gene_id)['gene_name']
            gene_mch_mmu_df[gene_label_mmu] = median_cluster_mch(get_gene_methylation('mouse_published', methylation_type, gene_id, level, True), level, cluster_type = 'cluster_ortholog')
        index = 0
        for gene_id in gene_id_list_hsa:
            if gene_id == None or gene_id == '':
                gene_label_hsa = "*N/A " + gene_mch_mmu_df.columns.values[index].upper() # ex. Cacna2d2 (mmu) -> CACNA2D2 (hsa)
                gene_mch_hsa_df[gene_label_hsa] = nan
            else:
                gene_label_hsa = get_gene_by_id('human_hv1_published', gene_id)['gene_name']
                gene_mch_hsa_df[gene_label_hsa] = median_cluster_mch(get_gene_methylation('human_hv1_published', methylation_type, gene_id, level, True), level, cluster_type = 'cluster_ortholog')
            index += 1
    else:
        gene_id_list_hsa = [ convert_gene_id_mmu_hsa(ensemble, gene_id) for gene_id in query.split() ]
        gene_id_list_mmu = [ find_orthologs(hsa_gene_id = gene_id)['mmu_gene_id'] for gene_id in gene_id_list_hsa ]
        for gene_id in gene_id_list_hsa:
            gene_label_hsa = get_gene_by_id('human_hv1_published', gene_id)['gene_name'] 
            gene_mch_hsa_df[gene_label_hsa] = median_cluster_mch(get_gene_methylation('human_hv1_published', methylation_type, gene_id, level, True), level, cluster_type = 'cluster_ortholog')
        index = 0
        for gene_id in gene_id_list_mmu:
            if gene_id == None or gene_id == '':
                gene_label_mmu = "*N/A " + gene_mch_hsa_df.columns.values[index].title() # ex. CACNA2D2 (hsa) -> Cacna2d2 (mmu)
                gene_mch_mmu_df[gene_label_mmu] = nan
            else:
                gene_label_mmu = get_gene_by_id('mouse_published', gene_id)['gene_name'] 
                gene_mch_mmu_df[gene_label_mmu] = median_cluster_mch(get_gene_methylation('mouse_published', methylation_type, gene_id, level, True), level, cluster_type = 'cluster_ortholog')
            index += 1

    if gene_mch_hsa_df.empty or gene_mch_mmu_df.empty:
        return FailToGraphException
    
    gene_mch_combined_df = gene_mch_hsa_df.join(gene_mch_mmu_df, how='outer')
    gene_mch_combined_df = gene_mch_combined_df[gene_mch_combined_df.index != '']

    if normalize_row:
        for gene in gene_mch_combined_df:
            gene_mch_combined_df[gene] = (gene_mch_combined_df[gene]-gene_mch_combined_df[gene].min()) / \
                    (gene_mch_combined_df[gene].max()-gene_mch_combined_df[gene].min())  

    if methylation_type == 'mch':
        methylation_type = 'mCH'
    else:
        methylation_type = 'mCG'

    title = "Orthologous gene body " + methylation_type + " by cluster"

    mch_mmu = [ gene_mch_combined_df[gene_name].tolist() for gene_name in gene_mch_mmu_df.columns.values ]
    mch_hsa = [ gene_mch_combined_df[gene_name].tolist() for gene_name in gene_mch_hsa_df.columns.values ]

    text_mmu, text_hsa, hover_mmu, hover_hsa = list(), list(), list(), list()
    for gene_name in gene_mch_mmu_df.columns.values:
        for cluster in gene_mch_combined_df.index:
            text_mmu.append(build_hover_text({
                    'Gene': gene_name,
                    'Cluster': cluster,
                    methylation_type: gene_mch_combined_df[gene_name][cluster],
                    })
                )   
        hover_mmu.append(text_mmu)
        text_mmu = []
    for gene_name in gene_mch_hsa_df.columns.values:
        for cluster in gene_mch_combined_df.index:
            text_hsa.append(build_hover_text({
                    'Gene': gene_name,
                    'Cluster': cluster,
                    methylation_type: gene_mch_combined_df[gene_name][cluster],
                    })
                )
        hover_hsa.append(text_hsa)
        text_hsa = []

    mch_combined = mch_mmu + mch_hsa
    flat_mch = list(chain.from_iterable(mch_combined))
    mch_dataframe = pd.DataFrame(flat_mch).dropna()
    start = mch_dataframe.quantile(0.05)[0].tolist()
    end = mch_dataframe.quantile(0.95).values[0].tolist()
    colorbar_tickval = list(arange(start, end, (end - start) / 4))
    colorbar_tickval[0] = start
    colorbar_tickval.append(end)
    colorbar_ticktext = [
        str(round(x, 3)) for x in arange(start, end, (end - start) / 4)
    ]
    if normalize_row == True:
        colorbar_ticktext[0] = str(round(start, 3))
    else:
        colorbar_ticktext[0] = '<' + str(round(start, 3))
    colorbar_ticktext.append('>' + str(round(end, 3)))

    # Due to a weird bug(?) in plotly, the number of elements in tickvals and ticktext 
    # must be greater than or equal to number of genes in query. Else, javascript throws 
    # Uncaught Typeerrors when trying to hover over genes. (Tomo 12/11/17)
    while len(colorbar_tickval) < len(gene_mch_hsa_df.columns):
        colorbar_tickval.insert(0,start)
        if normalize_row == True:
            colorbar_ticktext.insert(0, str(round(start, 3)))
        else:
            colorbar_ticktext.insert(0, '<' + str(round(start, 3)))

    trace_hsa = Heatmap(
            y=list(gene_mch_hsa_df.columns.values),
            x=gene_mch_combined_df.index,
            xaxis='x_hsa',
            yaxis='y_hsa',
            z=mch_hsa,
            text=hover_hsa,
            colorscale='Viridis',
            showscale=True,
            colorbar={
                'x': 1.0,
                'len': 1,
                'title': level.capitalize() + ' ' + methylation_type,
                'titleside': 'right',
                'tickmode': 'array',
                'tickvals': colorbar_tickval,
                'ticktext': colorbar_ticktext,
                'thickness': 10,
                'tickfont': {'size': 10}
                },
            hoverinfo='text'
            )
    trace_mmu = Heatmap(
            y=list(gene_mch_mmu_df.columns.values),    # Use hsa gene names to have same Y-axes for both
            x=gene_mch_combined_df.index,
            xaxis='x_mmu',
            yaxis='y_mmu',
            z=mch_mmu,
            text=hover_mmu,
            colorscale='Viridis',
            showscale=False,
            hoverinfo='text'
            )

    layout = Layout(
            title=title,
            titlefont={'color': 'rgba(1,2,2,1)',
                       'size': 16},
            autosize=True,
            height=450,

            hovermode='closest'
            )

    # Available colorscales:
    # https://community.plot.ly/t/what-colorscales-are-available-in-plotly-and-which-are-the-default/2079
    updatemenus = list([
        dict(
            buttons=list([
                dict(
                    args=['colorscale', 'Viridis'],
                    label='Viridis',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Bluered'],
                    label='Bluered',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Blackbody'],
                    label='Blackbody',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Electric'],
                    label='Electric',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Earth'],
                    label='Earth',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Jet'],
                    label='Jet',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Rainbow'],
                    label='Rainbow',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Picnic'],
                    label='Picnic',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'Portland'],
                    label='Portland',
                    method='restyle'
                ),
                dict(
                    args=['colorscale', 'YlGnBu'],
                    label='YlGnBu',
                    method='restyle'
                )
            ]),
            direction='down',
            showactive=True,
            x=0,
            xanchor='left',
            y=1.17,
            yanchor='top'
        )
    ])

    layout['updatemenus'] = updatemenus

    fig = tools.make_subplots(
            rows=1,
            cols=2,
            print_grid=False,
            shared_xaxes=True,
            shared_yaxes=False,
            subplot_titles=("Human_published", "Mouse_published"),
            )
    fig.append_trace(trace_hsa, row=1, col=1)
    fig.append_trace(trace_mmu, row=1, col=2)

    fig['layout'].update(layout)
    fig['layout']['xaxis1'].update({
        'side': 'bottom',
        'tickangle': -45,
        'tickfont': {
            'size': 12
             },
        'domain':[0, 0.45],
         })

    fig['layout']['xaxis2'].update({
        'side': 'bottom',
        'tickangle': -45,
        'tickfont': {
            'size': 12
             },
        'domain':[0.55, 1.0],
         })

    fig['layout']['yaxis1'].update({
        'visible':True,
        'tickangle': 15,
        'tickfont': {
            'size': 12
            },
        })

    fig['layout']['yaxis2'].update({
        'visible':True,
        'tickangle': 15,
        'tickfont': {
            'size': 12
            },
        })

    return plotly.offline.plot(
        figure_or_data=fig,
        output_type='div',
        show_link=False,
        include_plotlyjs=False)


@cache.memoize(timeout=3600)
def get_mch_box(ensemble, methylation_type, gene, grouping, clustering, level, outliers):
    """Generate gene body mCH box plot.

    Traces are grouped by cluster.

    Arguments:
        ensemble (str): Name of ensemble.
        methylation_type (str): Type of methylation to visualize.        "mch" or "mcg"
        gene (str):  Ensembl ID of gene for that ensemble.
        clustering (str): Different clustering algorithms and parameters. 'lv' = Louvain clustering.
        grouping (str): Variable to group cells by. "cluster", "annotation".
        level (str): "original" or "normalized" methylation values.
        outliers (bool): Whether if outliers should be displayed.

    Returns:
        str: HTML generated by Plot.ly.
    """
    gene = convert_gene_id_mmu_hsa(ensemble, gene)
    points = get_gene_methylation(ensemble, methylation_type, gene, grouping, clustering, level, outliers)
    context = methylation_type[1:]

    if points is None:
        raise FailToGraphException
    if points['annotation_'+clustering].nunique() <= 1:
        grouping = "cluster"
        print("**** Using cluster numbers")

    traces = OrderedDict()
    unique_groups = points[grouping+'_'+clustering].unique()
    max_cluster = len(unique_groups)

    colors = generate_cluster_colors(max_cluster, grouping)
    for point in points.to_dict('records'):
        name_prepend = ""
        if grouping == "cluster":
            name_prepend="cluster_"
        color = colors[int(np.where(unique_groups==point[grouping+'_'+clustering])[0]) % len(colors)]
        group = point[grouping+'_'+clustering]
        trace = traces.setdefault(group, Box(
                y=list(),
                name=name_prepend + str(group),
                marker={
                    'color': color,
                    'outliercolor': color,
                    'size': 6
                },
                boxpoints='suspectedoutliers',
                visible=True,
                showlegend=False,
                ))
        trace['y'].append(point[methylation_type + '/' + context + '_' + level])

    gene_name = get_gene_by_id(ensemble, gene)['gene_name']

    layout = Layout(
        autosize=True,
        height=450,
        title='Gene body ' + methylation_type + ' in each cluster: ' + gene_name,
        titlefont={'color': 'rgba(1,2,2,1)',
                   'size': 20},
#        legend={
#            'orientation': 'h',
#            'y': -0.3,
#            'traceorder': 'normal',
#        },
        xaxis={
            'title': 'Cluster',
            'titlefont': {
                'size': 17
            },
            'type': 'category',
            'anchor': 'y',
            'ticks': 'outside',
            'ticklen': 4,
            'tickangle': -45,
            'tickwidth': 0.5,
            'showticklabels': True,
            'tickfont': {
                'size': 12
            },
            'showline': True,
            'zeroline': False,
            'showgrid': True,
            'linewidth': 1,
            'mirror': True,
        },
        yaxis={
            'title': gene_name + ' ' + level.capitalize() + ' ' + methylation_type,
            'titlefont': {
                'size': 15
            },
            'type': 'linear',
            'anchor': 'x',
            'ticks': 'outside',
            # 'tickcolor': 'white',
            'ticklen': 4,
            'tickwidth': 0.5,
            'showticklabels': True,
            'tickfont': {
                'size': 12
            },
            'showline': True,
            'zeroline': False,
            'showgrid': True,
            'linewidth': 1,
            'mirror': True,
        },
    )

    return plotly.offline.plot(
        {
            'data': list(traces.values()),
            'layout': layout
        },
        output_type='div',
        show_link=False,
        include_plotlyjs=False)


@cache.memoize(timeout=3600)
def get_mch_box_two_ensemble(methylation_type, gene_mmu, gene_hsa, level, outliers):
    """Generate gene body mCH box plot for two ensemble.

    Traces are grouped by cluster and ordered by mm_hs_homologous_cluster.txt.
    Mouse clusters red, human clusters black.

    Arguments:
        methylation_type (str): Type of methylation to visualize.        "mch" or "mcg"
        gene_mmu (str):  Ensembl ID of gene mouse.
        gene_hsa (str):  Ensembl ID of gene human.
        level (str): "original" or "normalized" methylation values.
        outliers (bool): Whether if outliers should be displayed.

    Returns:
        str: HTML generated by Plot.ly.
    """
    gene_hsa = convert_gene_id_mmu_hsa('human_hv1_published', gene_hsa)
    gene_mmu = convert_gene_id_mmu_hsa('mouse_published', gene_mmu)
    points_mmu = get_gene_methylation('mouse_published', methylation_type, gene_mmu, level, outliers)
    points_hsa = get_gene_methylation('human_hv1_published', methylation_type, gene_hsa, level, outliers)
    cluster_order = get_ortholog_cluster_order()
    if points_mmu is None or points_hsa is None or not cluster_order is None:
        raise FailToGraphException

    gene_name = get_gene_by_id('mouse_published', gene_mmu)['gene_name']

    # EAM - This organizes the box plot into groups
    traces_mmu = Box(
        y=list(i.get(level) for i in points_mmu if i.get('cluster_ortholog')),
        x=list(i.get('cluster_ortholog') for i in points_mmu if i.get('cluster_ortholog')),
        marker={'color': 'red', 'outliercolor': 'red'},
        boxpoints='suspectedoutliers')
        
    traces_hsa = Box(
        y=list(i.get(level) for i in points_hsa if i.get('cluster_ortholog')),
        x=list(i.get('cluster_ortholog') for i in points_hsa if i.get('cluster_ortholog')),
        marker={'color': 'black', 'outliercolor': 'black'},
        boxpoints='suspectedoutliers')
    traces_combined = [traces_mmu, traces_hsa]

    layout = Layout(
        boxmode='group',
        autosize=True,
        height=450,
        showlegend=False,
        title='Gene body ' + methylation_type + ' in each cluster: ' + gene_name,
        titlefont={'color': 'rgba(1,2,2,1)',
                   'size': 20},
        # legend={
        #     'orientation': 'h',
        #     'x': -0.1,
        #     'y': -0.6,
        #     'traceorder': 'normal',
        # },
        xaxis={
            'title': '',
            'titlefont': {
                'size': 14
            },
            'type': 'category',
            'anchor': 'y',
            'ticks': 'outside',
            'tickcolor': 'rgba(51,51,51,1)',
            'ticklen': 4,
            'tickwidth': 0.5,
            'tickangle': -35,
            'showticklabels': True,
            'tickfont': {
                'size': 12
            },
            'showline': False,
            'zeroline': False,
            'showgrid': True,
        },
        yaxis={
            'title': gene_name+' '+level.capitalize() + ' mCH',
            'titlefont': {
                'size': 15
            },
            'type': 'linear',
            'anchor': 'x',
            'ticks': 'outside',
            'tickcolor': 'rgba(51,51,51,1)',
            'ticklen': 4,
            'tickwidth': 0.5,
            'showticklabels': True,
            'tickfont': {
                'size': 12
            },
            'showline': False,
            'zeroline': False,
            'showgrid': True,
        },
        shapes=[
            {
                'type': 'rect',
                'fillcolor': 'transparent',
                'line': {
                    'color': 'rgba(115, 115, 115, 1)',
                    'width': 1,
                    'dash': False
                },
                'yref': 'paper',
                'xref': 'paper',
                'x0': 0,
                'x1': 1,
                'y0': 0,
                'y1': 1
            },
        ],
        annotations=[{
            'text': '<b>■</b> Mouse',
            'x': 0.4,
            'y': 1.02,
            'ax': 0,
            'ay': 0,
            'showarrow': False,
            'font': {
                'color': 'red',
                'size': 12
            },
            'xref': 'paper',
            'yref': 'paper',
            'xanchor': 'left',
            'yanchor': 'bottom',
            'textangle': 0,
        }, {
            'text': '<b>■</b> Human',
            'x': 0.5,
            'y': 1.02,
            'ax': 0,
            'ay': 0,
            'showarrow': False,
            'font': {
                'color': 'Black',
                'size': 12
            },
            'xref': 'paper',
            'yref': 'paper',
            'xanchor': 'left',
            'yanchor': 'bottom',
            'textangle': 0,
        }])

    return plotly.offline.plot(
        {
            'data': traces_combined,
            'layout': layout
        },
        output_type='div',
        show_link=False,
        include_plotlyjs=False)
