"""Functions used to generate content. """
import datetime
import json
import math
import sys
import time

from collections import OrderedDict, Counter
from itertools import groupby, chain
from random import sample

import colorsys
import colorlover as cl
from flask import Blueprint, current_app, request
from sqlalchemy import exc, text
import numpy as np
from numpy import nan, linspace, arange, random
import pandas as pd
import plotly
from plotly import tools
from plotly.graph_objs import Layout, Annotation, Box, Scatter, Scatter3d, Heatmap, Bar # NOTE: Scattergl has some bugs
import sqlite3
from sqlite3 import Error
import plotly.figure_factory as ff
from scipy.spatial.distance import pdist, squareform
from scipy.cluster import hierarchy
from multiprocessing import Pool

from . import cache, db
from os import path

content = Blueprint('content', __name__) # Flask "bootstrap"

cluster_annotation_order = ['mL2/3', 'mL4', 'mL5-1', 'mL5-2', 'mDL-1', 'mDL-2', \
							'mL6-1', 'mL6-2', 'mDL-3', 'mVip', 'mNdnf-1', \
							'mNdnf-2', 'mPv', 'mSst-1', 'mSst-2', 'None']
methylation_types_order = ['mCH', 'mCG', 'mCA', 'mCHmCG', 'mCHmCA', 'mCAmCG']

num_sigfigs_ticklabels = 2;
ncells_max = 5000; # Max number of cells to show for scatter/box plots
log_file='/var/www/scmdb_py_dev/scmdb_log'

class FailToGraphException(Exception):
	"""Fail to generate data or graph due to an internal error."""
	pass

# @content.route('/content/metadata/')
# def get_metadata():
# 	result = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM cells LIMIT 1;").fetchall()
# 	result = [dict(r) for r in result]
# 	return json.dumps({"data": result})


@content.route('/content/ensembles')
def get_ensembles_summary():
	""" Retrieve data to be displayed in the "Ensembles" summary tabular page. 
		"/tabular/ensemble"
	"""
	regions = request.args.get('region', '').split()
	regions_lower = [ region.lower() for region in regions ]
	
	ensemble_list=[]
	# ensemble_list = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles").fetchall()
	ensemble_list = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles").fetchall()
	ensemble_list_atac = db.get_engine(current_app, 'snATAC_data').execute("SELECT * FROM ensembles").fetchall()
	# ensemble_list = ensemble_list.join(ensemble_list_atac, on="ensemble_id", how="outer")
	ensemble_list_ids = [ensemble['ensemble_id'] for ensemble in ensemble_list]
	for ensemble_atac in ensemble_list_atac:
		if (ensemble_atac['ensemble_id'] not in ensemble_list_ids):
			ensemble_list.append(ensemble_atac)

	total_methylation_cell_each_dataset = db.get_engine(current_app, 'methylation_data').execute("SELECT dataset, COUNT(*) as `num` FROM cells GROUP BY dataset").fetchall()
	total_methylation_cell_each_dataset = [ {d['dataset']: d['num']} for d in total_methylation_cell_each_dataset ]
	total_methylation_cell_each_dataset = { k.split('_',maxsplit=1)[1]: v for d in total_methylation_cell_each_dataset for k, v in d.items() }

	ensembles_cell_counts = []
	for ensemble in ensemble_list:
		ensemble_tbl = 'Ens' + str(ensemble['ensemble_id'])
		query_methylation = "SELECT dataset, COUNT(*) as `num` FROM cells INNER JOIN {} ON cells.cell_id = {}.cell_id GROUP BY dataset".format(ensemble_tbl, ensemble_tbl)
		try:
			methylation_cell_counts = db.get_engine(current_app, 'methylation_data').execute(query_methylation).fetchall()
			methylation_cell_counts = [ {d['dataset']: d['num']} for d in methylation_cell_counts]
			methylation_cell_counts = { k.split('_',maxsplit=1)[1]: v for d in methylation_cell_counts for k, v in d.items() }
		except:
			methylation_cell_counts = None

		query_snATAC = "SELECT dataset, COUNT(*) as `num` FROM cells INNER JOIN {} ON cells.cell_id = {}.cell_id GROUP BY dataset".format(ensemble_tbl, ensemble_tbl)
		try:
			snATAC_cell_counts = db.get_engine(current_app, 'snATAC_data').execute(query_snATAC).fetchall()
			snATAC_cell_counts = [ {d['dataset']: d['num']} for d in snATAC_cell_counts]
			snATAC_cell_counts = { k.split('_',maxsplit=1)[1]: v for d in snATAC_cell_counts for k, v in d.items() }
		except exc.ProgrammingError as e:
			snATAC_cell_counts = None

		annoj_exists = ensemble_annoj_exists(ensemble['ensemble_id'])
		ensembles_cell_counts.append( {"id": ensemble['ensemble_id'], 
									   "ensemble": ensemble['ensemble_name'], 
									   "ens_methylation_counts": methylation_cell_counts,
									   "ens_snATAC_counts": snATAC_cell_counts,
									   "public_access": ensemble['public_access'],
									   "description": ensemble['description'],
									   "annoj_exists": annoj_exists})

	ensembles_json_list = []
	for ens in ensembles_cell_counts:
		total_methylation_cells = 0
		total_snATAC_cells = 0
		datasets_in_ensemble_cell_count = []
		datasets_in_ensemble = []
		snATAC_datasets_in_ensemble = []
		ens_dict = {}
		if ens['ens_methylation_counts'] is not None:
			for dataset, count in ens['ens_methylation_counts'].items():
				total_methylation_cells += count
				datasets_in_ensemble.append('CEMBA_'+dataset)
				datasets_in_ensemble_cell_count.append(dataset+" ("+str(count)+" cells)")
				ens_dict[dataset] = str(count) + '/' + str(total_methylation_cell_each_dataset[dataset])
		if ens['ens_snATAC_counts'] is not None:
			for dataset, count in ens['ens_snATAC_counts'].items():
				total_snATAC_cells += count
				datasets_in_ensemble.append('CEMBA_'+dataset)
				snATAC_datasets_in_ensemble.append(dataset+" ("+str(count)+" cells)")

		if total_methylation_cells>200 or total_snATAC_cells>200: # Do not display ensembles that contain less than 200 total cells. (mainly RS2 data)

			ens_dict["ensemble_id"] = ens['id']
			ens_dict["ensemble_name"] = ens['ensemble']
			ens_dict["description"] = ens['description']
			ens_dict["datasets_rs1"] = ",  ".join(sorted([x for x in datasets_in_ensemble_cell_count if 'RS2' not in x]))
			ens_dict["datasets_rs2"] = ",  ".join(sorted([x for x in datasets_in_ensemble_cell_count if 'RS2' in x]))
			rs2_datasets_in_ensemble = sorted([x for x in datasets_in_ensemble if 'RS2' in x])
			ens_dict["target_regions_rs2_acronym"] = ""
			ens_dict["target_regions_rs2_descriptive"] = ""

			if len(rs2_datasets_in_ensemble) != 0:
				target_regions_query = "SELECT DISTINCT datasets.target_region, ABA_regions.ABA_description \
					FROM datasets \
					INNER JOIN ABA_regions ON ABA_regions.ABA_acronym=datasets.target_region \
					AND datasets.dataset in (" + ",".join(("%s",) * len(rs2_datasets_in_ensemble)) + ")"
				target_regions_result = db.get_engine(current_app, 'methylation_data').execute(target_regions_query, tuple(rs2_datasets_in_ensemble,)).fetchall()
				ens_dict["target_regions_rs2_acronym"] = ", ".join([ x.target_region for x in target_regions_result ])
				ens_dict["target_regions_rs2_descriptive"] = ", ".join([ x.ABA_description for x in target_regions_result ])

			ens_dict["snATAC_datasets_rs1"] = ",  ".join(sorted([x for x in snATAC_datasets_in_ensemble if 'RS2' not in x]))
			ens_dict["snATAC_datasets_rs2"] = ",  ".join(sorted([x for x in snATAC_datasets_in_ensemble if 'RS2' in x]))
			ens_dict["num_datasets"] = len(datasets_in_ensemble_cell_count)+len(snATAC_datasets_in_ensemble)

			slices_list_rs1 = [d.split('_')[0] for d in datasets_in_ensemble_cell_count if 'RS2' not in d]
			slices_list_rs1.extend([d.split('_')[0] for d in snATAC_datasets_in_ensemble if 'RS2' not in d])
			slices_list_rs2 = [d.split('_')[1][2:4] for d in datasets_in_ensemble_cell_count if 'RS2' in d]
			slices_list_rs2.extend([d.split('_')[0] for d in snATAC_datasets_in_ensemble if 'RS2' in d])
			slices_set = set(slices_list_rs1)
			slices_set.update(slices_list_rs2)
			ens_dict["slices"] = ",  ".join(sorted(list(slices_set)))
			ens_dict["total_methylation_cells"] = total_methylation_cells
			ens_dict["total_snATAC_cells"] = total_snATAC_cells

			if slices_set:
				ens_regions_query = "SELECT DISTINCT(ABA_acronym), ABA_description FROM ABA_regions WHERE `code` IN (" + ",".join("'"+x+"'" for x in list(slices_set)) + ")"
				ens_regions_result = db.get_engine(current_app, 'methylation_data').execute(ens_regions_query).fetchall()
				ens_regions_acronyms = [d['ABA_acronym'] for d in ens_regions_result]
				ens_regions_descriptions = [d['ABA_description'] for d in ens_regions_result]
				ens_dict["ABA_regions_acronym"] = ", ".join(ens_regions_acronyms).replace('+',', ')
				ens_dict["ABA_regions_description"] = ", ".join(ens_regions_descriptions).replace('+',', ')

			if ens['public_access'] == 0:
				ens_dict["public_access_icon"] = "fas fa-lock"
				ens_dict["public_access_color"] = "black"
			else:
				ens_dict["public_access_icon"] = "fas fa-lock-open"
				ens_dict["public_access_color"] = "green"

			ens_dict["annoj_exists"] = ens['annoj_exists']


			if regions == ['None']:
				ensembles_json_list.append(ens_dict)
			else:
				for region in regions_lower:
					if region in ens_dict["ABA_regions_acronym"].lower():
						ensembles_json_list.append(ens_dict)
						break

	ens_json = json.dumps(ensembles_json_list)

	return ens_json


@content.route('/content/datasets/<rs>')
def get_datasets_summary(rs):
	""" Retrieve data to be displayed in the RS1 and RS2 summmary tabular page. 
		"/tabular/dataset/rs1"
		"/tabular/dataset/rs2"
	
		Arguments:
			rs = Research Segment. Either "rs1" or "rs2"
	"""
	
	if rs == "rs1":
		dataset_list = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM datasets WHERE dataset NOT LIKE 'CEMBA_RS2_%%'").fetchall()
		dataset_list += db.get_engine(current_app, 'snATAC_data').execute("SELECT * FROM datasets WHERE dataset NOT LIKE 'CEMBA_RS2_%%'").fetchall()
		
		# This is a hack to get unique values in a list of dictionaries
		dataset_list = list({x['dataset']:x for x in dataset_list}.values()); 
		total_methylation_cell_each_dataset = db.get_engine(current_app, 'methylation_data').execute("SELECT dataset, COUNT(*) as `num` FROM cells WHERE dataset NOT LIKE 'CEMBA_RS2_%%' GROUP BY dataset").fetchall()
		total_snATAC_cell_each_dataset = db.get_engine(current_app, 'snATAC_data').execute("SELECT dataset, COUNT(*) as `num` FROM cells WHERE dataset NOT LIKE 'CEMBA_RS2_%%' GROUP BY dataset").fetchall()
	elif rs == "rs2":
		dataset_list = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM datasets WHERE dataset LIKE 'CEMBA_RS2_%%'").fetchall()
		total_methylation_cell_each_dataset = db.get_engine(current_app, 'methylation_data').execute("SELECT dataset, COUNT(*) as `num` FROM cells WHERE dataset LIKE 'CEMBA_RS2_%%' GROUP BY dataset").fetchall()
		total_snATAC_cell_each_dataset = db.get_engine(current_app, 'snATAC_data').execute("SELECT dataset, COUNT(*) as `num` FROM cells WHERE dataset LIKE 'CEMBA_RS2_%%' GROUP BY dataset").fetchall()
	elif rs == "all":
		dataset_list = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM datasets").fetchall()
		dataset_list += db.get_engine(current_app, 'snATAC_data').execute("SELECT * FROM datasets").fetchall()
		# This is a hack to get unique values in a list of dictionaries
		dataset_list = list({x['dataset']:x for x in dataset_list}.values()); 
		total_methylation_cell_each_dataset = db.get_engine(current_app, 'methylation_data').execute("SELECT dataset, COUNT(*) as `num` FROM cells GROUP BY dataset").fetchall()
		total_snATAC_cell_each_dataset = db.get_engine(current_app, 'snATAC_data').execute("SELECT dataset, COUNT(*) as `num` FROM cells GROUP BY dataset").fetchall()
	else:
		return

	total_methylation_cell_each_dataset = [ {d['dataset']: d['num']} for d in total_methylation_cell_each_dataset ]
	total_methylation_cell_each_dataset = { k: v for d in total_methylation_cell_each_dataset for k, v in d.items() }
	total_snATAC_cell_each_dataset = [ {d['dataset']: d['num']} for d in total_snATAC_cell_each_dataset ]
	total_snATAC_cell_each_dataset = { k: v for d in total_snATAC_cell_each_dataset for k, v in d.items() }

	dataset_cell_counts = []
	for dataset in dataset_list:
		try:
			num_snATAC_cells = total_snATAC_cell_each_dataset[dataset['dataset']]
		except KeyError as e:
			num_snATAC_cells = 0

		if "RS2" not in dataset['dataset']:
			brain_region_code = dataset['dataset'].split('_')[1]
			research_segment = "RS1"
		else:
			brain_region_code = dataset['dataset'].split('_')[2]
			brain_region_code = brain_region_code[-2:]
			research_segment = "RS2"
			
		regions_sql = db.get_engine(current_app, 'methylation_data').execute("SELECT ABA_description FROM ABA_regions WHERE ABA_acronym=%s", (dataset['brain_region'],)).fetchone()
		if regions_sql is not None:
			ABA_regions_descriptive = regions_sql['ABA_description'].replace('+', ', ')
		else: 
			ABA_regions_descriptive = ""

		if rs == "rs1":
			try:
				dataset_cell_counts.append( {"dataset_name": dataset['dataset'], 
											 "sex": dataset['sex'],
											 "methylation_cell_count": total_methylation_cell_each_dataset[dataset['dataset']],
											 "snATAC_cell_count": num_snATAC_cells,
											 "ABA_regions_acronym": dataset['brain_region'].replace('+', ', '),
											 "ABA_regions_descriptive": ABA_regions_descriptive,
											 "slice": brain_region_code,
											 "date_added": str(dataset['date_online']),
											 "description": dataset['description'] })
			except:
				dataset_cell_counts.append( {"dataset_name": dataset['dataset'], 
											 "sex": dataset['sex'],
											 "methylation_cell_count": 0,
											 "snATAC_cell_count": num_snATAC_cells,
											 "ABA_regions_acronym": dataset['brain_region'].replace('+', ', '),
											 "ABA_regions_descriptive": ABA_regions_descriptive,
											 "slice": brain_region_code,
											 "date_added": str(dataset['date_online']),
											 "description": dataset['description'] })
		else:
			target_region_sql = db.get_engine(current_app, 'methylation_data').execute("SELECT ABA_description FROM ABA_regions WHERE ABA_acronym=%s", (dataset['target_region'],)).fetchone()
			if target_region_sql is not None:
				target_region_descriptive = target_region_sql['ABA_description'].replace('+', ', ')
			else:
				target_region_descriptive = ""

			try:
				dataset_cell_counts.append( {"dataset_name": dataset['dataset'],
											 "research_segment": research_segment,
											 "sex": dataset['sex'],
											 "methylation_cell_count": total_methylation_cell_each_dataset[dataset['dataset']],
											 "snATAC_cell_count": num_snATAC_cells,
											 "ABA_regions_acronym": dataset['brain_region'].replace('+', ', '),
											 "ABA_regions_descriptive": ABA_regions_descriptive,
											 "slice": brain_region_code,
											 "date_added": str(dataset['date_online']),
											 "description": dataset['description'],
											 "target_region_acronym": dataset['target_region'],
											 "target_region_descriptive": target_region_descriptive})
			except:
				dataset_cell_counts.append( {"dataset_name": dataset['dataset'],
											 "research_segment": research_segment,
											 "sex": dataset['sex'],
											 "methylation_cell_count": 0,
											 "snATAC_cell_count": num_snATAC_cells,
											 "ABA_regions_acronym": dataset['brain_region'].replace('+', ', '),
											 "ABA_regions_descriptive": ABA_regions_descriptive,
											 "slice": brain_region_code,
											 "date_added": str(dataset['date_online']),
											 "description": dataset['description'],
											 "target_region_acronym": dataset['target_region'],
											 "target_region_descriptive": target_region_descriptive})

	return json.dumps(dataset_cell_counts)

@content.route("/content/check_ensembles/<new_ensemble_name>/<new_ensemble_datasets>")
def check_ensemble_similarities(new_ensemble_name, new_ensemble_datasets):
	"""
	Used by the "request_new_ensemble" page. Checks if the new ensemble has any similarities with pre-existing ensembles to prevent duplication of ensembles.
	"""

	existing_ensembles = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles").fetchall()
	existing_ensembles_list = [ dict(d) for d in existing_ensembles ]
	existing_ensembles_names_list = [ d['ensemble_name'] for d in existing_ensembles ]

	# New ensemble_name must be unique.
	if new_ensemble_name in existing_ensembles_names_list:
		return json.dumps({"result": "failure", "reason": "The name {} is already in use, please choose a different name".format(new_ensemble_name)})

	new_ensemble_datasets = new_ensemble_datasets.split('+')
	
	query = "SELECT cell_id FROM cells WHERE dataset IN (" + ",".join(('%s',)*len(new_ensemble_datasets)) + ")"
	cells_in_new_ensemble = db.get_engine(current_app, 'methylation_data').execute(query, tuple(new_ensemble_datasets)).fetchall()
	cells_in_new_ensemble_set = set([ cell['cell_id'] for cell in cells_in_new_ensemble ])

	if len(cells_in_new_ensemble_set) <= 200:
		return json.dumps({"result": "failure", "reason": "Ensembles must contain more than 200 cells."})
	
	same_datasets_in_both = []
	new_ensemble_datasets_set = set(new_ensemble_datasets)
	for existing_ensemble in existing_ensembles_list:
		existing_ensemble_datasets = set(existing_ensemble['datasets'].split(','))
		datasets_difference = new_ensemble_datasets_set ^ existing_ensemble_datasets # datasets in new or existing but not both (difference).
		if len(datasets_difference) == 0:
			same_datasets_in_both.append(existing_ensemble)

	for similar_ensemble in same_datasets_in_both:
		query = "SELECT cell_id FROM Ens{}".format(similar_ensemble['ensemble_id'])
		cells_in_similar_ensemble = db.get_engine(current_app, 'methylation_data').execute(query).fetchall()
		cells_in_similar_ensemble_set = set([ cell['cell_id'] for cell in cells_in_similar_ensemble ])
		different_cells = cells_in_new_ensemble_set ^ cells_in_similar_ensemble_set

		# If a pre-existing ensemble with the same datasets also has the same exact cells as the new ensemble, tell user a duplicate ensemble exists
		if len(different_cells) == 0:
			return json.dumps({"result": "failure", "reason": "Another ensemble with the same cells already exists: {}.".format(similar_ensemble['ensemble_name'])})

	# If none of the pre-existing ensembles with the same datasets has the same exact cells as the new ensemble, warn user that similar ensembles exist.
	if len(same_datasets_in_both) > 0:
		return json.dumps({"result": "warning", "reason": "The following pre-existing ensembles are similar: "+", ".join(("%s",)*len(same_datasets_in_both)) %(tuple([ ensemble['ensemble_name'] for ensemble in same_datasets_in_both]))+". Are you sure you want to request the new ensemble?"})

	# Success
	return json.dumps({"result": "success", 
					   "reason": "Click submit to finalize request.",
					   "new_ensemble_name": new_ensemble_name, 
					   "new_ensemble_datasets": new_ensemble_datasets, 
					   "new_ensemble_cells": list(cells_in_new_ensemble_set)})

# Utilities
@cache.memoize(timeout=1800)
def ensemble_exists(ensemble, modality='methylation'):
	"""Check if data for a given ensemble exists 

	Arguments:
		ensemble (str): Name of ensemble.
		modality (str): methylation, snATAC or RNA

	Returns:
		bool: Whether if given ensemble exists
	"""

	if modality=='methylation':
		ensemble_name = 'ensemble_id'
	else:
		ensemble_name = 'snmc_ensemble_id'
	result = db.get_engine(current_app, modality+'_data').execute("SELECT * FROM ensembles WHERE "+ensemble_name+"=%s", (str(ensemble),)).fetchone()

	if result is None:
		return 0
	else:
		return 1

# Utilities
@cache.memoize(timeout=1800)
def ensemble_annoj_exists(ensemble):
	"""Check if AnnoJ browser for a given ensemble exists 

	Arguments:
		ensemble (str): Name of ensemble.
		
	Returns:
		bool: Whether if given ensemble exists
	"""

	ensemble = str(ensemble)
	ensemble = ensemble.strip('Ens')
	result = path.isfile('/var/www/html/annoj_private/CEMBA/browser/fetchers/mc_cemba/mc_single_merged_mCG_cluster_mCHmCG_lv_npc50_k30_1_Ens'+str(ensemble)+'.php');

	return result

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
		if v is not None:
			text += '{k}: {v}<br>'.format(k=k, v=str(v))

	return text.strip('<br>')


def generate_cluster_colors(num, grouping):
	"""Generate a list of colors given number needed.

	Arguments:
		num (int): Number of colors needed. n <= 35.

	Returns:
		list: strings containing RGB-style strings e.g. rgb(255,255,255).
	"""

	if (grouping == 'dataset' or grouping == 'target_region') and num > 2 and num <= 9:
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
def get_genes_of_module(module):
	"""Generates list of genes in selected module.
	Arguments:
		module (str): Name of module to query for.
	Returns:
		Dataframe of gene_name and gene_id of each gene in the module for the corresponding
	"""

	modules_result = db.get_engine(current_app, 'methylation_data').execute("SELECT module, mmu_gene_id, mmu_gene_name FROM gene_modules WHERE module=%s", (module,)).fetchall()
	genes_in_module = [ {'module': d['module'], 'gene_id': d['mmu_gene_id'], 'gene_name': d['mmu_gene_name']} for d in modules_result ]

	return genes_in_module

@cache.memoize(timeout=3600)
def get_cluster_marker_genes(ensemble, clustering):
	"""Retrieves list of top marker genes for each cluster in a clustering of an ensemble.
	Arguments:
		ensemble (str): Ensemble id. ie. Ens0, Ens1...
		clustering (str): Different clustering algorithms and parameters. 'lv' = Louvain clustering. ie. mCH_lv_npc50_k5
	Returns:
	list of dicts. ie [{'clustering': 'mCH_lv_npc50_k5', 'cluster': 1, 'rank': 1, 'gene_id': 'ENSMUSG_########', 'gene_name': 'Gad2'}]
	"""

	if ';' in ensemble:
		return []

	query = "SELECT clustering, cluster, rank, genes.gene_id, genes.gene_name \
		FROM {0}_cluster_marker_genes \
		INNER JOIN genes ON {0}_cluster_marker_genes.gene_id = genes.gene_id \
		WHERE clustering = %s".format(ensemble)

	try:
		result = db.get_engine(current_app, 'methylation_data').execute(query, (clustering,)).fetchall()
	except exc.ProgrammingError as e:
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_cluster_marker_genes): {}".format(str(now), e))
		sys.stdout.flush()
		return []

	result = [ dict(x) for x in result ]
	num_genes = result[-1]['rank']
	num_clusters = result[-1]['cluster']
	columns = [ 'cluster_'+str(i+1) for i in range(num_clusters) ] 

	rows = []
	for i in range(num_genes):
		row = {}
		row = dict(zip(columns, [ gene['gene_name'] for gene in result[i:: num_genes] ] ))
		row['rank'] = i+1
		rows.append(row)


	to_json = {'columns': columns}
	to_json['rows'] = rows

	return to_json

@cache.memoize(timeout=3600)
def median_cluster_mch(gene_info, grouping, clustering):
	"""Returns median mch level of a gene for each cluster.

		Arguments:
			gene_info (dict): mCH data for each sample. Keys are samp(cell), tsne_x, tsne_y, cluster_label, cluster_ordered, original, normalized.
			level (str): "original" or "normalized" methylation values.

		Returns:
			dict: Cluster_label (key) : median mCH level (value).
	"""

	if gene_info is not None:
		if grouping == 'annotation':
			gene_info.fillna({'annotation_'+clustering: 'None'}, inplace=True)
			return gene_info.groupby('annotation_'+clustering, sort=False)[gene_info.columns[-1]].median()
		elif grouping == 'cluster':
			gene_info.fillna({'cluster_'+clustering: 'None'}, inplace=True)
			return gene_info.groupby('cluster_'+clustering, sort=False)[gene_info.columns[-1]].median()
		elif grouping == 'dataset':
			gene_info.fillna({'dataset': 'None'}, inplace=True)
			return gene_info.groupby('dataset', sort=False)[gene_info.columns[-1]].median()
		elif grouping == 'target_region':
			gene_info['target_region'].fillna('N/A', inplace=True)
			return gene_info.groupby('target_region', sort=False)[gene_info.columns[-1]].median()
		elif grouping == 'slice':
			datasets_all_cells = gene_info['dataset'].tolist()
			slices_list = [d.split('_')[1] if 'RS2' not in d else d.split('_')[2][2:4] for d in datasets_all_cells]
			gene_info['slice'] = slices_list
			return gene_info.groupby('slice', sort=False)[gene_info.columns[-2]].median()
		elif grouping == 'sex':
			return gene_info.groupby('sex', sort=False)[gene_info.columns[-1]].median()
		else:
			return None
	else:
		return None

@cache.memoize(timeout=3600)
def mean_cluster(gene_info, grouping, modality='ATAC'):
	"""Returns median mch level of a gene for each cluster.

		Arguments:
			gene_info (dict): mCH data for each sample. Keys are samp(cell), tsne_x, tsne_y, cluster_label, cluster_ordered, original, normalized.
			level (str): "original" or "normalized" methylation values.
			modality (str): 'ATAC','RNA'
		Returns:
			dict: Cluster_label (key) : mean mCH level (value).
	"""

	if grouping == 'annotation':
		gene_info.fillna({'annotation_'+modality: 'None'}, inplace=True)
	if grouping != 'dataset':
		return gene_info.groupby(grouping+'_'+modality, sort=False)['normalized_counts'].mean()
	else:
		return gene_info.groupby(grouping, sort=False)['normalized_counts'].mean()

	if grouping == 'annotation':
		gene_info.fillna({'annotation_'+modality: 'None'}, inplace=True)
		return gene_info.groupby('annotation_'+modality, sort=False)['normalized_counts'].mean()
	elif grouping == 'cluster':
		return gene_info.groupby('cluster_'+modality, sort=False)['normalized_counts'].mean()
	elif grouping == 'dataset':
		return gene_info.groupby('dataset', sort=False)['normalized_counts'].mean()
	elif grouping == 'target_region':
		gene_info['target_region'].fillna('N/A', inplace=True)
		return gene_info.groupby('target_region', sort=False)['normalized_counts'].mean()
	else:
		return None

@cache.memoize(timeout=3600)
def get_ensemble_info(ensemble_name=str(), ensemble_id=str()):
	"""
	Gets information regarding an ensemble. Requires either the ensemble name or ensemble id
	"""
	
	if ensemble_name:
		# result = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles WHERE ensemble_name=%s", (ensemble_name,)).fetchone()
		result = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles WHERE ensemble_name=%s", (ensemble_name,)).fetchone()
		if result is None:
			result = db.get_engine(current_app, 'snATAC_data').execute("SELECT * FROM ensembles WHERE ensemble_name=%s", (ensemble_name,)).fetchone()
	else:
		ensemble_id = int(filter(str.isdigit, ensemble_id))
		result = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM ensembles WHERE ensemble_id=%s", (ensemble_id,)).fetchone()
		if result is None:
			result = db.get_engine(current_app, 'snATAC_data').execute("SELECT * FROM ensembles WHERE ensemble_name=%s", (ensemble_id,)).fetchone()	

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
		print("[{}] ERROR in app(get_methylation_tsne_options): {}".format(str(now), e))
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

	list_algorithms_clustering = list(set([x.split('_')[1] for x in list_clustering_types]))

	#list_mc_types_clustering = sorted(list(set([x.split('_')[0] for x in list_clustering_types])), key=lambda mC_type: methylation_types_order.index(mC_type))
	#list_algorithms_clustering = sorted(list(set([x.split('_')[1] for x in list_clustering_types if list_mc_types_clustering[0] == x.split('_')[0]])))
	#list_npc_clustering = sorted(list(set([int(x.split('_')[2].replace('npc', '')) for x in list_clustering_types if (list_mc_types_clustering[0]+'_'+list_algorithms_clustering[0]) in x])))
	#list_k_clustering = sorted(list(set([int(x.split('_')[3].replace('k', '')) for x in list_clustering_types if (list_mc_types_clustering[0]+'_'+list_algorithms_clustering[0]+'_npc'+str(list_npc_clustering[0])) in x])))

	query = "SELECT * FROM cells LIMIT 1"
	df_metadata = pd.read_sql(query, con=db.get_engine(current_app, 'methylation_data'))
	df_metadata = df_metadata.drop(list(df_metadata.filter(regex='cell_.*', axis='columns')), axis=1)
	list_metadata = ['cluster','annotation']
	list_metadata += list([i for i in df_metadata.columns.values])

	return {'all_tsne_settings': list_tsne_types, 
			'tsne_methylation': list_mc_types_tsne,
			'all_clustering_settings': list_clustering_types,
			'all_clustering_settings2': dict_clustering_types_and_numclusters,
			'clustering_algorithms': list_algorithms_clustering,
			'tsne_dimensions': list_dims_tsne_first,
			'tsne_perplexity': list_perp_tsne_first,
			'methylation_metadata_fields': list_metadata,}

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
		print("[{}] ERROR in app(get_snATAC_tsne_options): {}".format(str(now), e))
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
def get_gene_by_name(gene_query):
	"""Retrieve gene information by name. Mainly used to fill gene search bar.
	Does not search for exact matches.

	Arguments:
		gene_query (list): List of gene name strings

	Returns:
		DataFrame: Info for queried gene(s). Columns are gene_id, gene_name, chr, start, end, strand, gene_type.
	"""

	gene_query = [ gene.lower()+"%" for gene in gene_query ]

	sql_query = "SELECT * FROM genes WHERE " + "lower(gene_name) LIKE %s OR " * len(gene_query)
	sql_query = sql_query[:-3]

	df = pd.read_sql(sql_query, params=(gene_query,), con=db.get_engine(current_app, 'methylation_data'))

	return df.to_dict('records')

@cache.memoize(timeout=3600)
def get_gene_by_name_exact(gene_query):
	"""Same as get_gene_by_name but for exact matches only.

	Arguments:
		gene_query (list): List of gene name strings

	Returns:
		DataFrame: Info for queried gene(s). Columns are gene_id, gene_name, chr, start, end, strand, gene_type.
	"""

	gene_query = [ gene.lower() for gene in gene_query ]
	placeholders_str = "%s, " * len(gene_query)
	placeholders_str = placeholders_str[:-2]
	sql_query = "SELECT * FROM genes WHERE " + "lower(gene_name) IN (" + placeholders_str+ ")"
	sql_query += " ORDER BY CASE lower(gene_name) "

	for i, gene in enumerate(gene_query):
		sql_query += "WHEN '{}' THEN {} ".format(gene, i+1)
	sql_query += "END"

	df = pd.read_sql(sql_query, params=(gene_query,), con=db.get_engine(current_app, 'methylation_data'))

	return df.to_dict('records')
	
@cache.memoize(timeout=3600)
def get_gene_by_id(gene_query):
	"""Retrieve gene information by gene id.

	Arguments:
		gene_query (list): list of gene_id strings.

	Returns:
		DataFrame: Info for queried gene. Columns are gene_id, gene_name, chr, start, end, strand, gene_type.
	"""

	gene_query_wildcard = [ gene+'%' for gene in gene_query ] 
	sql_query = "SELECT * FROM genes WHERE " + "gene_id LIKE %s OR " * len(gene_query_wildcard)
	sql_query = sql_query[:-3]

	df = pd.read_sql(sql_query, params=(gene_query_wildcard,), con=db.get_engine(current_app, 'methylation_data'))

	#reorder genes to original order since SQL doesn't keep order.
	new_index = []
	for index, row in df.iterrows():
		for i, gene_id in enumerate(gene_query):
			if gene_id in row['gene_id']:
				new_index.append(i)
				break
	df.index = new_index
	df.sort_index(inplace=True)

	return df.to_dict('records')

@cache.memoize(timeout=3600)
def get_corr_genes(ensemble, query):
	"""Get correlated genes of a certain gene of a ensemble. 
	
		Arguments:
			ensemble(str): Ensemble identifier. (Eg. Ens0, Ens1, Ens2...).
			query(str): Gene ID.
		
		Returns:
			dict: information of genes that are correlated with target gene.
	"""
	if ";" in query:
		return []

	try:
		corr_genes = db.get_engine(current_app, 'methylation_data').execute("SELECT * FROM {}_correlated_genes WHERE gene1 LIKE %s".format(ensemble), (query+'%%',)).fetchall()
	except exc.ProgrammingError as e:
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_corr_genes): {}".format(str(now), e))
		sys.stdout.flush()
		return []

	corr_genes = [ {"rank": i+1, "gene_name": get_gene_by_id(row.gene2)[0]['gene_name'], "correlation": row.correlation, "gene_id": row.gene2} for i, row in enumerate(corr_genes)]
	return corr_genes

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
	if grouping in ['annotation','cluster']:
		groupingu = ensemble+"."+grouping+"_"+clustering
	elif grouping in ['NeuN']:
		groupingu = "CONCAT('NeuN',cells."+grouping+")"
	else:
		groupingu = "cells."+grouping

	if 'ndim2' in tsne_type:
		query = "SELECT cells.cell_id, cells.dataset, %(ensemble)s.cluster_%(clustering)s, datasets.target_region, \
			%(ensemble)s.annotation_%(clustering)s, %(gene_table_name)s.%(methylation_type)s, \
			cells.global_%(methylation_type)s, %(groupingu)s as grouping, \
			%(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, \
			%(gene_table_name)s.%(context)s, datasets.sex \
			FROM cells \
			INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
			LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
			LEFT JOIN datasets ON cells.dataset = datasets.dataset" % {'ensemble': ensemble, 'groupingu': groupingu,
																	   'gene_table_name': gene_table_name,
																	   'tsne_type': tsne_type,
																	   'methylation_type': methylation_type,
																	   'context': context,
																	   'clustering': clustering,}
		# query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, cells.global_%(methylation_type)s, \
		# 	%(ensemble)s.annotation_%(clustering)s, %(ensemble)s.cluster_%(clustering)s, \
		# 	%(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, \
		# 	%(gene_table_name)s.%(methylation_type)s, %(gene_table_name)s.%(context)s, \
		# 	datasets.target_region, datasets.sex \
		# 	FROM cells \
		# 	INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
		# 	LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
		# 	LEFT JOIN datasets ON cells.dataset = datasets.dataset" % {'ensemble': ensemble,
		# 															   'gene_table_name': gene_table_name,
		# 															   'tsne_type': tsne_type,
		# 															   'methylation_type': methylation_type,
		# 															   'context': context,
		# 															   'clustering': clustering,}
	else:
		query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, cells.global_%(methylation_type)s, \
			%(ensemble)s.annotation_%(clustering)s, %(ensemble)s.cluster_%(clustering)s, \
			%(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, %(ensemble)s.tsne_z_%(tsne_type)s, \
			%(gene_table_name)s.%(methylation_type)s, %(gene_table_name)s.%(context)s, \
			datasets.target_region, datasets.sex \
			FROM cells \
			INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
			LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
			LEFT JOIN datasets ON cells.dataset = datasets.dataset" % {'ensemble': ensemble, 
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
		df.fillna({'grouping': 'None'}, inplace=True)
		df['annotation_cat'] = pd.Categorical(df['grouping'], cluster_annotation_order)
		df.sort_values(by='annotation_cat', inplace=True)
		df.drop('annotation_cat', axis=1, inplace=True)
	elif grouping == 'cluster':
		df.sort_values(by='cluster_'+clustering, inplace=True)

	# # Add number of cells in each cluster to the cluster name
	# ncells = df.groupby(by=grouping, sort=False).count()
	# for label, i in df[grouping]:
	#     ncells[0]

	return df

def get_gene_from_mysql(ensemble, gene_table_name, methylation_type, clustering, tsne_type, grouping='cluster'):
	"""Helper function to fetch a gene's methylation information from mysql.

	TODO: Don't need to fetch tsne info, annotations etc. except once

	Returns:
		DataFrame
	"""

	context = methylation_type[1:]
	if grouping in ['annotation','cluster']:
		groupingu = ensemble+"."+grouping+"_"+clustering
	else:
		groupingu = "cells."+grouping

	t0=datetime.datetime.now()
	# print(' Running get_gene_from_mysql for '+gene_table_name+' : '+str(t0)+'; ', file=open(log_file,'a'))# EAM - Profiling SQL
	if tsne_type=='noTSNE':
		query = "SELECT %(gene_table_name)s.%(methylation_type)s, %(gene_table_name)s.%(context)s, \
			FROM %(ensemble)s  \
			LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble, 
																	   'gene_table_name': gene_table_name,
																	   'methylation_type': methylation_type,
																	   'context': context,}
	elif 'ndim2' in tsne_type:
		query = "SELECT cells.cell_id, cells.dataset, %(ensemble)s.cluster_%(clustering)s, datasets.target_region, \
			%(ensemble)s.annotation_%(clustering)s, %(gene_table_name)s.%(methylation_type)s, \
			cells.global_%(methylation_type)s, %(groupingu)s as grouping, \
			%(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, \
			%(gene_table_name)s.%(context)s, datasets.sex \
			FROM cells \
			INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
			LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
			LEFT JOIN datasets ON cells.dataset = datasets.dataset" % {'ensemble': ensemble, 'groupingu': groupingu,
																	   'gene_table_name': gene_table_name,
																	   'tsne_type': tsne_type,
																	   'methylation_type': methylation_type,
																	   'context': context,
																	   'clustering': clustering,}
		
		# query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, cells.global_%(methylation_type)s, \
		# 	%(ensemble)s.annotation_%(clustering)s, %(ensemble)s.cluster_%(clustering)s, \
		# 	%(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, \
		# 	%(gene_table_name)s.%(methylation_type)s, %(gene_table_name)s.%(context)s, \
		# 	datasets.target_region, datasets.sex \
		# 	FROM cells \
		# 	INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
		# 	LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
		# 	LEFT JOIN datasets ON cells.dataset = datasets.dataset" % {'ensemble': ensemble, 
		# 															   'gene_table_name': gene_table_name,
		# 															   'tsne_type': tsne_type,
		# 															   'methylation_type': methylation_type,
		# 															   'context': context,
		# 															   'clustering': clustering,}
	else: # 3D tSNE
		query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, cells.global_%(methylation_type)s, \
			%(ensemble)s.annotation_%(clustering)s, %(ensemble)s.cluster_%(clustering)s, \
			%(ensemble)s.tsne_x_%(tsne_type)s, %(ensemble)s.tsne_y_%(tsne_type)s, %(ensemble)s.tsne_z_%(tsne_type)s, \
			%(gene_table_name)s.%(methylation_type)s, %(gene_table_name)s.%(context)s, \
			datasets.target_region, datasets.sex \
			FROM cells \
			INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
			LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
			LEFT JOIN datasets ON cells.dataset = datasets.dataset" % {'ensemble': ensemble, 
																	   'gene_table_name': gene_table_name,
																	   'tsne_type': tsne_type,
																	   'methylation_type': methylation_type,
																	   'context': context,
																	   'clustering': clustering,}
	try:
		df = pd.read_sql(query, db.get_engine(current_app, 'methylation_data'))
	except exc.ProgrammingError as e:
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_mult_gene_methylation): {}".format(str(now), e))
		sys.stdout.flush()
		return None

	return df


@cache.memoize(timeout=3600)
def get_mult_gene_methylation(ensemble, methylation_type, genes, grouping, clustering, level, tsne_type):
	"""Return averaged methylation data ponts for a set of genes.

	Data from ID-to-Name mapping and tSNE points are combined for plot generation.

	Arguments:
		ensemble (str): Name of ensemble.
		methylation_type (str): Type of methylation to visualize. "mch", "mcg", etc
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

	# This query is just to fix gene id's missing the Ensembl version number. 
	# Necessary because the table name must match exactly with whats on the MySQL database.
	# Ex. ENSMUSG00000026787 is fixed to ENSMUSG00000026787.3
	first_query = "SELECT gene_id FROM genes WHERE gene_id LIKE %s" + " OR gene_id LIKE %s" * (len(genes)-1)
	result = db.get_engine(current_app, 'methylation_data').execute(first_query, (genes,)).fetchall()

	gene_table_names = ['gene_' + gene_id[0].replace('.','_') for gene_id in result]

	df_all = pd.DataFrame()
	
	############
	# t0=datetime.datetime.now()
	# print(' Starting mysql queries '+str(t0)+'; ', file=open(log_file,'a'))# EAM - Profiling SQL
	# print('Pool size 12', file=open(log_file,'a'))
	# with Pool(12) as pool:
	# 	df_list = [pool.apply_async(get_gene_from_mysql,
 #                                       args=(ensemble, gene_table_name, methylation_type, clustering, tsne_type)).get()
 #                               		for gene_table_name in gene_table_names ]
	# t1=datetime.datetime.now()
	# print(' Finished mysql queries '+str(t1)+'; ', file=open(log_file,'a'))# EAM - Profiling SQL
	# print(' Time '+str(t1-t0)+'; ', file=open(log_file,'a'))# EAM - Profiling SQL

	t0=datetime.datetime.now()
#	print(' Starting mysql queries '+str(t0)+'; ', file=open(log_file,'a'))# EAM - Profiling SQL
#	print('Pool size 1', file=open(log_file,'a'))
	df_list = []
	for i, gene_table_name in enumerate(gene_table_names):
		t0a=datetime.datetime.now()
		if i>0:
			tsne_typeu='noTSNE'
		else:
			tsne_typeu=tsne_type
		df_all = df_all.append(get_gene_from_mysql(ensemble, gene_table_name, methylation_type, clustering, tsne_typeu, grouping))
		if i==0:
			df_coords=df_all
		t1a=datetime.datetime.now()
		# print(str(i)+' : ',str(t1a-t0a), file=open(log_file,'a'))


	t1=datetime.datetime.now()
#	print(' Finished mysql queries '+str(t1)+'; ', file=open(log_file,'a'))# EAM - Profiling SQL
#	print(' Time '+str(t1-t0)+'; ', file=open(log_file,'a'))# EAM - Profiling SQL
	############

	df_all[[methylation_type, context]] = df_all[[methylation_type, context]].apply(pd.to_numeric)

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
		df_coords.sort_values(by='annotation_cat', inplace=True)
		df_coords.drop('annotation_cat', axis=1, inplace=True)
	elif grouping == 'cluster':
		df_coords.sort_values(by='cluster_'+clustering, inplace=True)

	return df_coords

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

	genes = genes_query.split()

	gene_name_str = ""
	x, y, text, mch = list(), list(), list(), list()

	if len(genes) == 1:
		points = get_gene_methylation(ensemble, methylation_type, genes[0], grouping, clustering, level, True, tsne_type)
		gene_name_str = get_gene_by_id([ genes[0] ])[0]['gene_name']
		title = 'Gene body ' + methylation_type + ': ' + gene_name_str
	else:
		points = get_mult_gene_methylation(ensemble, methylation_type, genes, grouping, clustering, level, tsne_type)
		gene_infos = get_gene_by_id(genes)
		for i, gene in enumerate(gene_infos):
			if i > 0 and i % 10 == 0:
				gene_name_str += '<br>'
			gene_name_str += gene['gene_name'] + '+'
		gene_name_str = gene_name_str[:-1]
		title = 'Avg. Gene body ' + methylation_type + ': <br>' + gene_name_str

	if points is None:
		raise FailToGraphException

	### TSNE ### 
	if grouping == 'annotation':
		if grouping+'_'+clustering not in points.columns or points[grouping+'_'+clustering].nunique() <= 1:
			grouping = "cluster"
			clustering = "mCH_lv_npc50_k30"
			# print("**** Using cluster_mCH_lv_npc50_k30")

	datasets = points['dataset'].unique().tolist()
	if grouping in ['cluster','annotation','dataset','NeuN','sex']:
		unique_groups = points['grouping'].unique().tolist()
	else:
		# For continuous (numerical) metadata (like global_mCH), don't use discrete clusters
		unique_groups = ['All cells']

	if grouping == 'cluster':
		annotation_additional_y = 0.025 # Necessary because legend items overlap with legend title (annotation) when there are many legend items
		num_clusters = points['cluster_'+clustering].max()
	else:
		annotation_additional_y = 0.00 
		num_clusters = len(unique_groups)
	
	colors = generate_cluster_colors(len(unique_groups), grouping)
	symbols = ['circle', 'square', 'cross', 'triangle-up', 'triangle-down', 'octagon', 'star', 'diamond']
	
	traces_tsne = OrderedDict()

	legend_x = -.17
	layout_width = 1100
	grouping_clustering = grouping
	if grouping == 'cluster' or grouping == 'annotation':
		layout_width = 1000
		legend_x = -.14
		grouping_clustering = grouping+'_'+clustering

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

	context = methylation_type[1:]

	## 2D tSNE coordinates ##
	if 'ndim2' in tsne_type:
		for i, group in enumerate(unique_groups):
			if group == 'All cells':
				# Continuous variable
				points_group = points
				color_num = i
				color = points['grouping']
			else:
				# Categorial variable
				points_group = points[points['grouping']==group]
				color_num = i
				color = colors[color_num]

			if grouping_clustering.startswith('cluster'):
				group_str = 'cluster_' + str(group)
			elif grouping_clustering== "dataset":
				group = "_".join(group.split('_')[1:])
				group_str = group
			else:
				group_str = group
			if group == 'All cells':
				group_str = ''
			
			
			trace2d = traces_tsne.setdefault(color_num, Scatter(
				x=list(),
				y=list(),
				text=list(),
				mode='markers',
				visible=True,
				name=group_str,
				legendgroup=group,
				marker={
					   'color': color,
					   'colorscale': 'Viridis',
					   'size': marker_size,
					   #'symbol': symbols[datasets.index(dataset)],
				},
				hoverinfo='text'))
			trace2d['x'] = points_group['tsne_x_'+tsne_type].values.tolist()
			trace2d['y'] = points_group['tsne_y_'+tsne_type].values.tolist()
			trace2d['text'] = [build_hover_text(OrderedDict([('Annotation', point[4]),
														  ('Cluster', point[2]),
														  ('RS2 Target Region', point[3]),
														  ('Dataset', point[1]),
														  ('<b>'+grouping+'</b>', point[7]),]))
							   for point in points_group.itertuples(index=False)]

		### METHYLATION SCATTER ### 
		x = points['tsne_x_' + tsne_type].tolist()
		y = points['tsne_y_' + tsne_type].tolist()
		mch = points[methylation_type + '/' + context + '_' + level]
		text_methylation = [build_hover_text(OrderedDict([('Annotation', point[4]),
														  ('Cluster', point[2]),
														  ('RS2 Target Region', point[3]),
														  ('Dataset', point[1]),
														  ('<b>'+level.title()+' '+methylation_type+'</b>', round(point[5], 6)),]))
							for point in points.itertuples(index=False)]


		mch_dataframe = pd.DataFrame(mch)
		start = mch_dataframe.dropna().quantile(ptile_start)[0].tolist()
		end = mch_dataframe.dropna().quantile(ptile_end).values[0].tolist()
		end = max(end,start+0.01)
		mch_colors = [set_color_by_percentile(x, start, end) for x in mch]

		colorbar_tickval = list(arange(start, end, (end - start) / 4))
		colorbar_tickval[0] = start
		colorbar_tickval.append(end)
		colorbar_ticktext = [
			str(round(x, num_sigfigs_ticklabels)) for x in arange(start, end, (end - start) / 4)
		]
		colorbar_ticktext[0] = '<' + str(round(start, num_sigfigs_ticklabels))
		colorbar_ticktext.append('>' + str(round(end, num_sigfigs_ticklabels)))

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
			height=550,
			width=layout_width,
			#title=title,
			#titlefont={'color': 'rgba(1,2,2,1)',
			#           'size': 16},
			legend={'x':legend_x,
					'y':0.95,
					'tracegroupgap': 0.5},
			margin={'l': 0,
					'r': 0,
					'b': 30,
					't': 130,},
			xaxis={
				'domain': [0, 0.49],
				'type': 'linear',
				'ticks': '',
				'dtick': 10,
				'tickwidth': 0,
				'showticklabels': False,
				'showline': True,
				'showgrid': False,
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
				'showline': True,
				'showgrid': False,
				'zeroline': False,
				'linecolor': 'black',
				'linewidth': 0.5,
				'mirror': False,
				'scaleanchor': 'y',
				'range':[bottom_x,top_x]
			},
			yaxis={
				'domain': [0, 1],
				'type': 'linear',
				'ticks': '',
				'dtick': 10,
				'tickwidth': 0,
				'showticklabels': False,
				'showline': True,
				'showgrid': False,
				'side': 'right',
				'zeroline': False,
				'linecolor': 'black',
				'linewidth': 0.5,
				'mirror': False,
				'range':[bottom_y,top_y]
			},
			hovermode='closest',)
			#hoverdistance='10',)

		
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
		fig['layout']['annotations'].extend([Annotation(text=grouping.title(),
														x=legend_x+0.05,
														y=1.02 + annotation_additional_y,
														xanchor="left",
														yanchor="top",
														showarrow=False,
														xref="paper",
														yref="paper",
														font={'size': 12,
															  'color': 'gray',})])

		fig['layout']['annotations'].extend([Annotation(text=title,
														x=0.5,
														y=1.3,
														xanchor="center",
														yanchor="top",
														showarrow=False,
														xref="paper",
														yref="paper",
														font={'size': 16,
															  'color': 'black',})])


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
			trace3d['text'] = [build_hover_text(OrderedDict([('Dataset', point[2]),
															 ('Annotation', point[4]),
															 ('Cluster', point[5]),]))
							   for point in points_group.itertuples(index=False)]

		### METHYLATION SCATTER ### 
		x = points['tsne_x_' + tsne_type].tolist()
		y = points['tsne_y_' + tsne_type].tolist()
		z = points['tsne_z_' + tsne_type].tolist()
		mch = points[methylation_type + '/' + context + '_' + level]
		text_methylation = [build_hover_text(OrderedDict([('Annotation', point[4]),
														  ('Cluster', point[5]),
														  ('<b>'+methylation_type+'</b>', round(point[-1], 6)),]))
							for point in points.itertuples(index=False)]


		mch_dataframe = pd.DataFrame(mch)
		start = mch_dataframe.dropna().quantile(ptile_start)[0].tolist()
		end = mch_dataframe.dropna().quantile(ptile_end).values[0].tolist()
		end = max(end,start+0.01)
		mch_colors = [set_color_by_percentile(x, start, end) for x in mch]

		colorbar_tickval = list(arange(start, end, (end - start) / 4))
		colorbar_tickval[0] = start
		colorbar_tickval.append(end)
		colorbar_ticktext = [
			str(round(x, num_sigfigs_ticklabels)) for x in arange(start, end, (end - start) / 4)
		]
		colorbar_ticktext[0] = '<' + str(round(start, num_sigfigs_ticklabels))
		colorbar_ticktext.append('>' + str(round(end, num_sigfigs_ticklabels)))

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
			width=1000,
			title=title,
			titlefont={'color': 'rgba(1,2,2,1)',
					   'size': 16},
			legend={'x':-.1, 'y':1},
			margin={'l': 49,
					'r': 0,
					'b': 30,
					't': 100,
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
	points = get_gene_methylation(ensemble, methylation_type, gene, grouping, clustering, level, outliers)
	context = methylation_type[1:]

	if points is None:
		raise FailToGraphException
	if grouping == 'annotation' and points['annotation_'+clustering].nunique() <= 1:
		grouping = "cluster"
		print("**** Using cluster numbers")

	traces = OrderedDict()
	if grouping == "dataset":
		unique_groups = points["dataset"].unique()
	elif grouping == 'target_region':
		points['target_region'].fillna('N/A', inplace=True)
		unique_groups = points['target_region'].unique()
	elif grouping == 'slice':
		datasets_all_cells = points['dataset'].tolist()
		slices_list = [d.split('_')[1] if 'RS2' not in d else d.split('_')[2][2:4] for d in datasets_all_cells]
		points['slice'] = slices_list
		slices_set = set(slices_list)
		unique_groups = np.array(list(slices_set))
	elif grouping == 'sex':
		unique_groups = points['sex'].unique()
	elif grouping == 'cluster' or grouping == 'annotation':
		unique_groups = points[grouping+'_'+clustering].unique()
	else:
		grouping = 'cluster'
		unique_groups = points[grouping+'_'+clustering].unique()
		# raise FailToGraphException
	num_clusters = len(unique_groups)

	colors = generate_cluster_colors(num_clusters, grouping)
	for point in points.to_dict('records'):
		name_prepend = ""
		if grouping == "dataset" or grouping == 'target_region' or grouping == 'slice' or grouping == 'sex':
			color = colors[int(np.where(unique_groups==point[grouping])[0]) % len(colors)]
			group = point[grouping]
		else:
			if grouping == "cluster":
				name_prepend="cluster_"
			color = colors[int(np.where(unique_groups==point[grouping+'_'+clustering])[0]) % len(colors)]
			group = point[grouping+'_'+clustering]
		if outliers:
			boxpoints='suspectedoutliers';
		else:
			boxpoints=False

		trace = traces.setdefault(group, Box(
				y=list(),
				name=name_prepend + str(group),
				marker={
					'color': color,
					'outliercolor': color,
					'size': 6
				},
				boxpoints=boxpoints,
				visible=True,
				showlegend=False,
				))
		trace['y'].append(point[methylation_type + '/' + context + '_' + level])

	gene_name = get_gene_by_id([ gene ])[0]['gene_name']

	layout = Layout(
		autosize=True,
		height=450,
		width=1000,
		title='Gene body ' + methylation_type + ' in each cluster: ' + gene_name,
		titlefont={'color': 'rgba(1,2,2,1)',
				   'size': 20},
		xaxis={
			'title': grouping.title(),
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
		normal_or_original = '(normalized by row)'
	else:
		normal_or_original = ''

	title = level.title() + " gene body " + methylation_type + " by cluster " + normal_or_original + ": <br>"
	genes = query.split()

	s=''
	for i in genes:
		s=s+','+i

	gene_labels = list()
	gene_info_df = pd.DataFrame()
	gene_infos = get_gene_by_id(genes)
	for i, gene in enumerate(gene_infos):
		gene_name = gene['gene_name']
		gene_labels.append(gene_name)
		if i > 0 and i % 10 == 0:
			title += "<br>"
		title += gene_name + "+"
		gene_info_df[gene_name] = median_cluster_mch(get_gene_methylation(ensemble, methylation_type, gene['gene_id'], grouping, clustering, level, True), grouping, clustering)
		if gene_info_df[gene_name].empty:
			raise FailToGraphException

	title = title[:-1] # Gets rid of last '+'

	gene_info_df.reset_index(inplace=True)
	if grouping == 'annotation':
		gene_info_df['annotation_cat'] = pd.Categorical(gene_info_df['annotation_'+clustering], cluster_annotation_order)
		gene_info_df.sort_values(by='annotation_cat', inplace=True)
		gene_info_df.drop('annotation_cat', axis=1, inplace=True)
		gene_info_df.set_index(grouping+'_'+clustering, inplace=True)
	elif grouping == 'cluster':
		gene_info_df.sort_values(by='cluster_'+clustering, inplace=True)
		gene_info_df.set_index(grouping+'_'+clustering, inplace=True)
	elif grouping == 'dataset' or grouping == 'target_region' or grouping == 'slice' or grouping == 'sex':
		gene_info_df.sort_values(by=grouping, inplace=True)
		gene_info_df.set_index(grouping, inplace=True)
	else:
		raise FailToGraphException

	# For some reason, Plotly doesn't allow 'None' as a group on the x-axis for heatmaps.
	if gene_info_df.index.tolist() == ['None']: 
		gene_info_df.index = ['N/A']

	clusters_labels = gene_info_df.index.tolist()
	if grouping == 'cluster':
		clusters_labels = ['Cluster '+str(i) for i in clusters_labels]

	normal_or_original = 'Original'
	if normalize_row:
		for gene in gene_info_df:
			# z-score
			# gene_info_df[gene] = (gene_info_df[gene] - gene_info_df[gene].mean()) / gene_info_df[gene].std()
			# min-max
			gene_range = gene_info_df[gene].max() - gene_info_df[gene].min()
			if (gene_range==0):
				gene_range = 1
			gene_info_df[gene] = (gene_info_df[gene] - gene_info_df[gene].min()) / gene_range
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
			text.append(build_hover_text(OrderedDict([('Gene', key),
													  (grouping.title(), x[j]),
													  (methylation_type, mch[i][j])
													  ])))
			j += 1
		hover.append(text)
		text = []
		i += 1

	flat_mch = list(chain.from_iterable(mch))
	mch_dataframe = pd.DataFrame(flat_mch).dropna()

	# Hierarchical clustering and dendrogram
	mch = np.array(mch)
	figure = ff.create_dendrogram(mch, orientation="right", labels=tuple([i for i in range(len(genes))]), 
		colorscale=['bbbbbbb']) # TODO: Figure out how to set the colorscale
	for i in range(len(figure['data'])):
		figure['data'][i]['xaxis'] = 'x2'
	dendro_leaves = figure['layout']['yaxis']['ticktext']
	dendro_leaves = list(map(int, dendro_leaves))
	mch = mch[dendro_leaves,:] # Reorder the genes according to the clustering
	genes_labels = [gene_labels[i] for i in dendro_leaves]
	hover_old = hover
	# hover = [hover_old[i] for i in dendro_leaves]
	hover = [str(i) for i in dendro_leaves]

	dendro_top = ff.create_dendrogram(mch.transpose(), orientation="bottom", labels=tuple([i for i in range(mch.shape[1])]), 
		colorscale=['bbbbbbb'])
	for i in range(len(dendro_top['data'])):
		dendro_top['data'][i]['yaxis'] = 'y2'
	dendro_top_leaves = dendro_top['layout']['xaxis']['ticktext']
	dendro_top_leaves = list(map(int, dendro_top_leaves))
	mch = mch[:,dendro_top_leaves] # Reorder the genes according to the clustering
	clusters_labels = [clusters_labels[i] for i in dendro_top_leaves]
	mch = list(mch)
	figure['data'].extend(dendro_top['data'])

	# Set color scale limits
	start = mch_dataframe.quantile(ptile_start).values[0].tolist()
	end = mch_dataframe.quantile(ptile_end).values[0].tolist()
	end = max(end,start+0.01)
	
	colorbar_tickval = list(arange(start, end, (end - start) / 4))
	colorbar_tickval[0] = start
	colorbar_tickval.append(end)
	colorbar_ticktext = [
		str(round(x, num_sigfigs_ticklabels)) for x in arange(start, end, (end - start) / 4)
	]
	if normalize_row == True:
		colorbar_ticktext[0] = str(round(start, num_sigfigs_ticklabels))
	else:
		if (round(start,num_sigfigs_ticklabels)) == 0:
			colorbar_ticktext[0] = str(round(start,num_sigfigs_ticklabels))
		else:
			colorbar_ticktext[0] = '<' + str(round(start, num_sigfigs_ticklabels))
	colorbar_ticktext.append('>' + str(round(end, num_sigfigs_ticklabels)))

	# Due to a weird bug(?) in plotly, the number of elements in tickvals and ticktext 
	# must be greater than or equal to number of genes in query. Else, javascript throws 
	# Uncaught Typeerrors when trying to hover over genes. (Tomo 12/11/17)
	while len(colorbar_tickval) < len(genes):
		colorbar_tickval.insert(0,start)
		if normalize_row == True:
			colorbar_ticktext.insert(0, str(round(start, num_sigfigs_ticklabels)))
		else:
			colorbar_ticktext.insert(0, '<' + str(round(start, num_sigfigs_ticklabels)))

	trace = Heatmap(
		x=dendro_top_leaves,
		y=dendro_leaves,
		z=mch,
		xtype="array", ytype="array",
		text=hover,
		colorscale='Viridis',
		colorbar={
			'x': 1.0,
			'len': 0.5,
			'title': level.capitalize() + ' ' + methylation_type,
			'titleside': 'right',
			'tickmode': 'array',
			'tickvals': colorbar_tickval,
			'ticktext': colorbar_ticktext,
			'thickness': 10,
			'tickfont': {'size': 10}
			},
		hoverinfo='text',
		zmin=start,zmax=end,zauto=False, # Clip the extreme edges of the colorscale
		)
	trace['y'] = figure['layout']['yaxis']['tickvals']
	trace['x'] = dendro_top['layout']['xaxis']['tickvals']
	figure['data'].extend([trace])

	layout = Layout(
		height=max(600*len(genes)/20,550), # EAM Adjust the height of the heatmap according to the number of genes displayed
		width=1000,
		paper_bgcolor='rgba(0,0,0,0)',
		plot_bgcolor='rgba(0,0,0,0)',
		showlegend=False,
		hovermode='closest',
		title=title,
		# titlefont={'color': 'rgba(1,2,2,1)',
		#            'size': 16},
		margin={'l': 0,
				'r': 0,
				'b': 100,
				't': 150,},
		xaxis={
			'side': 'bottom',
			'tickangle': -45,
			'title': 'Clusters',
			'tickfont': {'size': 12},
			'showticklabels': True,
			'tickmode': 'array',
			'tickvals':trace['x'],
			'ticktext':clusters_labels,
			},
		yaxis={
			# 'tickangle': 15,
			'tickfont': {'size': 12},
			'showticklabels': True,
			'ticks':"outside",
			'tickmode': 'array',
			'tickvals':trace['y'],
			'ticktext':genes_labels,
			},
		)
	layout['yaxis'].update({'domain': [0, .85]})
	layout['xaxis'].update({'domain': [0.2, 1]})
	layout.update({'hovermode': 'closest'})
	layout.update({'xaxis2': {
			'showticklabels': False
			}})
	layout.update({'yaxis2': {
			'showticklabels': False
			}})
	layout['xaxis2'].update({'domain': [0, 0.1]})
	layout['yaxis2'].update({'domain': [0.86, 1]})
	for xx in ['xaxis','yaxis','xaxis2','yaxis2']:
		layout[xx].update({'mirror': False,
						   'showgrid': False,
						   'showline': False,
						   'zeroline': False})

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
			x=-0.1,
			xanchor='left',
			y=1.43,
			yanchor='top'
		)
	])

	layout['updatemenus'] = updatemenus

	# layout['annotations'].extend([Annotation(text=title,
	#                                          x=0.5,
	#                                          y=1.3,
	#                                          xanchor="center",
	#                                          yanchor="top",
	#                                          showarrow=False,
	#                                          xref="paper",
	#                                          yref="paper",
	#                                          font={'size': 16,
	#                                                'color': 'black',})])

	figure['layout'] = layout

	return plotly.offline.plot(figure,
		output_type='div',
		show_link=False,
		include_plotlyjs=False)

@cache.memoize(timeout=3600)
def get_clusters(ensemble, grouping, clustering):
	"""Return information about all the clusters

	Arguments:
		ensemble (str): Name of ensemble.
		grouping (str): Variable for grouping cells. "cluster", "annotation", or "dataset".
		clustering (str): Different clustering algorithms and parameters. 'lv' = Louvain clustering.

	Returns:
		DataFrame
	"""

	# Prevent SQL injected since column names cannot be parameterized.
	if ";" in ensemble or ";" in grouping or ";" in clustering:
		return None

	if grouping in ['annotation','cluster']:
		groupingu = ensemble+"."+grouping+"_"+clustering
	elif grouping in ['NeuN']:
		groupingu = "CONCAT('NeuN',cells."+grouping+")"
	else:
		groupingu = "cells."+grouping

	# Get methylation info
	query = "SELECT count(cells.cell_id) ncells, 'snmC' as modality, \
		%(groupingu)s as groups \
		FROM cells \
		INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
		GROUP BY groups " % {'ensemble': ensemble,
					'groupingu': groupingu,
					'clustering': clustering}
	try:
		df = pd.read_sql(query, db.get_engine(current_app, 'methylation_data'))
	except exc.ProgrammingError as e:
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_clusters): {}".format(str(now), e))
		sys.stdout.flush()
		return None

	# Get snATAC info
	query = "SELECT count(cells.cell_id) ncells, 'snATAC' AS modality, %(ensemble)s.cluster_ATAC groups \
		FROM cells \
		INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
		GROUP BY groups " % {'ensemble': ensemble,
					'grouping': grouping,
					'clustering': clustering}

	try:
		df_atac = pd.read_sql(query, db.get_engine(current_app, 'snATAC_data'))
		df=df.append(df_atac)
	except exc.ProgrammingError as e:
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_clusters): {}".format(str(now), e))
		sys.stdout.flush()


	# Get snRNA info
	query = "SELECT count(cells.cell_id) ncells, 'RNA' AS modality, %(ensemble)s.cluster_RNA groups \
		FROM cells \
		INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
		GROUP BY groups " % {'ensemble': ensemble,
					'grouping': grouping,
					'clustering': clustering}

	try:
		df_rna = pd.read_sql(query, db.get_engine(current_app, 'RNA_data'))
		df=df.append(df_rna)
	except exc.ProgrammingError as e:
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_clusters): {}".format(str(now), e))
		sys.stdout.flush()

	return df

@cache.memoize(timeout=3600)
def get_clusters_bar(ensemble, grouping, clustering, normalize, outliers):
	"""Generate clusters bar plot.

	Traces are grouped by modality (mch, ATAC).

	Arguments:
		ensemble (str): Name of ensemble.
		clustering (str): Different clustering algorithms and parameters. 'lv' = Louvain clustering.
		grouping (str): Variable to group cells by. "cluster", "annotation".
		outliers (bool): Whether if outliers should be displayed.

	Returns:
		str: HTML generated by Plot.ly.
	"""
	if grouping not in ['cluster','annotation','dataset','NeuN']:
		grouping = 'cluster'
	clusters = get_clusters(ensemble, grouping, clustering)

	if (normalize=='true'):
		clusters['ncells_norm'] = clusters.groupby('groups')['ncells'].transform(lambda x: 100*x / x.sum())
		clusters['y'] = clusters['ncells_norm']
		ytitle = 'Percent of cells per cluster'
	else:
		clusters['y'] = clusters['ncells']
		ytitle = 'Number of cells per cluster'

	# Stacked bar chart by modality. Appropriate for integrated clustering only
	mu = clusters['modality'].unique().tolist()
	data = list();
	for mi in mu:
		clustersu = clusters[clusters['modality']==mi]
		trace = Bar(
			y=clustersu['y'],
			x=clustersu['groups'],
			name=mi+' cells',
			hoverinfo='text',
			)
		if (normalize=='true'):
			trace['text'] = [str(round(i,1))+'% '+mi+' cells' for i in clustersu['y']]
		else:
			trace['text'] = [str(i)+' '+mi+' cells' for i in clustersu['y']]
		data.append(trace)
	
	layout = Layout(
	    autosize=True,
	    height=450,
	    width=1000,
	    title=ytitle,
	    titlefont={'color': 'rgba(1,2,2,1)',
	               'size': 20},
	    barmode='stack',
	    xaxis={
	        'title': 'Cluster',
	        'titlefont': {
	            'size': 17
	        },
	        'type': 'category',
	        'tickvals':clusters['groups'].unique(),
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
	        'title': "Number of cells",
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
			'data': data,
			'layout': layout
		},
		output_type='div',
		show_link=False,
		include_plotlyjs=False)

### TODO: Refactor the code to combine the ATAC and RNA into one set of functions...
### snATAC
@cache.memoize(timeout=3600)
def get_gene_snATAC(ensemble, gene, grouping, outliers, smoothing=False):
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
	
	if smoothing:
		counts_type='smoothed_normalized_counts'
	else:
		counts_type='normalized_counts'

	query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, \
		%(ensemble)s.annotation_ATAC, %(ensemble)s.cluster_ATAC, \
		%(ensemble)s.tsne_x_ATAC, %(ensemble)s.tsne_y_ATAC, \
		%(gene_table_name)s.%(counts_type)s as normalized_counts, \
		datasets.target_region \
		FROM cells \
		INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
		LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
		LEFT JOIN datasets ON cells.dataset = datasets.dataset \
		ORDER BY RAND() LIMIT %(ncells)s" % {'ensemble': ensemble, 
																   'gene_table_name': gene_table_name,
																   'counts_type': counts_type,
																   'ncells': ncells_max}

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
		df.sort_values(by='annotation_cat', inplace=True)
		df.drop('annotation_cat', axis=1, inplace=True)
	elif grouping == 'cluster':
		df.sort_values(by='cluster_ATAC', inplace=True)

	df['normalized_counts'].fillna(0, inplace=True)
	
	return df

def get_gene_snatac_from_mysql(ensemble, gene_table_name, counts_type, tsne_type):
	"""Helper function to fetch a gene's snatac information from mysql.

	Returns:
		DataFrame
	"""

	t0=datetime.datetime.now()
	if tsne_type=='noTSNE':
		query = "SELECT %(gene_table_name)s.%(counts_type)s as normalized_counts \
			FROM %(ensemble)s  \
			LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
			LIMIT 5000" % {'ensemble': ensemble, 
			   'gene_table_name': gene_table_name,
			   'counts_type': counts_type,}
	else:
		query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, \
			%(ensemble)s.annotation_ATAC, %(ensemble)s.cluster_ATAC, \
			%(ensemble)s.tsne_x_ATAC, %(ensemble)s.tsne_y_ATAC, \
			%(gene_table_name)s.%(counts_type)s as normalized_counts, \
			datasets.target_region \
			FROM cells \
			INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
			LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
			LEFT JOIN datasets ON cells.dataset = datasets.dataset \
			LIMIT 5000" % {'ensemble': ensemble, 
																	   'gene_table_name': gene_table_name,
																	   'counts_type': counts_type}
	try:
		df = pd.read_sql(query, db.get_engine(current_app, 'snATAC_data'))
	except exc.ProgrammingError as e:
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_gene_snatac_from_mysql): {}".format(str(now), e))
		sys.stdout.flush()
		return None

	t1=datetime.datetime.now()
#	print(' Running get_gene_snatac_from_mysql for '+gene_table_name+' : '+str(t1-t0)+'; ', file=open(log_file,'a')) # EAM - Profiling SQL
#	print(' query: '+query, file=open(log_file,'a'))

	return df

@cache.memoize(timeout=1800)
def get_mult_gene_snATAC(ensemble, genes, grouping, smoothing=False):
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
	
	if smoothing:
		counts_type='smoothed_normalized_counts'
	else:
		counts_type='normalized_counts'

	# first = True
	# for i,gene_table_name in enumerate(gene_table_names):
	# 	t0=datetime.datetime.now()
	# 	if first:
	# 		query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, \
	# 			%(ensemble)s.annotation_ATAC, %(ensemble)s.cluster_ATAC, \
	# 			%(ensemble)s.tsne_x_ATAC, %(ensemble)s.tsne_y_ATAC, \
	# 			%(gene_table_name)s.%(counts_type)s as normalized_counts, \
	# 			datasets.target_region \
	# 			FROM cells \
	# 			INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
	# 			LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
	# 			LEFT JOIN datasets ON cells.dataset = datasets.dataset" % {'ensemble': ensemble, 
	# 																	   'gene_table_name': gene_table_name,
	# 																	   'counts_type': counts_type}
	# 	else:
	# 		query = "SELECT %(gene_table_name)s.%(counts_type)s as normalized_counts \
	# 			FROM %(ensemble)s  \
	# 			LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble, 
	# 																	   'gene_table_name': gene_table_name,
	# 																	   'counts_type': counts_type}

	# 		# query = "SELECT %(gene_table_name)s.%(counts_type)s as normalized_counts \
	# 		# 	FROM %(gene_table_name)s  \
	# 		# 	INNER JOIN %(ensemble)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id" % {'ensemble': ensemble, 
	# 		# 															   'gene_table_name': gene_table_name,
	# 		# 															   'counts_type': counts_type}			

	# 	try:
	# 		df_all = df_all.append(pd.read_sql(query, db.get_engine(current_app, 'snATAC_data')))
	# 	except exc.ProgrammingError as e:
	# 		now = datetime.datetime.now()
	# 		print("[{}] ERROR in app(get_mult_gene_snATAC): {}".format(str(now), e))
	# 		sys.stdout.flush()
	# 		return None
		
	# 	if first:
	# 		df_coords = df_all
	# 	first = False
	# 	t1=datetime.datetime.now()
	# 	print(str(i)+' loading snATAC: '+str(t1-t0), file=open(log_file, 'a'))

	t0=datetime.datetime.now()
#	print('Pool size 1', file=open(log_file,'a'))
	df_list = []
	df = get_gene_snatac_from_mysql(ensemble, gene_table_names[0], counts_type, 'TSNE')
	df_all = df_all.append(df)
	df_coords = df;
	# df_list.append()
	for gene_table_name in gene_table_names[1:]:
		df_all = df_all.append(get_gene_snatac_from_mysql(ensemble, gene_table_name, counts_type, 'noTSNE'))
	# with Pool(12) as pool:
	# 	df_list = [pool.apply_async(get_gene_snatac_from_mysql,
 # 	                                      args=(ensemble, gene_table_name, counts_type, 'noTSNE')).get()
 #                               		for gene_table_name in gene_table_names[1:] ]

	t1=datetime.datetime.now()
#	print('All done: '+str(t1-t0), file=open(log_file,'a'))
	# df_coords=df_list[0]
	# for df in df_list:
		# df_all = df_all.append(df)

	t1=datetime.datetime.now()
#	print('All done: '+str(t1-t0), file=open(log_file,'a'))

	if df_all.empty: # If no data in column, return None 
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_gene_snATAC): No snATAC data for {}".format(str(now), ensemble))
		sys.stdout.flush()
		return None

	df_all['normalized_counts'].fillna(0, inplace=True)

	df_avg_methylation = df_all.groupby(by='cell_id', as_index=False)['normalized_counts'].mean()
	df_coords.update(df_avg_methylation)

	if grouping == 'annotation':
		df_coords.fillna({'annotation_ATAC': 'None'}, inplace=True)
		df_coords['annotation_cat'] = pd.Categorical(df_coords['annotation_ATAC'], cluster_annotation_order)
		df_coords.sort_values(by='annotation_cat', inplace=True)
		df_coords.drop('annotation_cat', axis=1, inplace=True)
	elif grouping == 'cluster':
		df_coords.sort_values(by='cluster_ATAC', inplace=True)
	return df_coords

@cache.memoize(timeout=1800)
def get_snATAC_scatter(ensemble, genes_query, grouping, ptile_start, ptile_end, tsne_outlier_bool, smoothing=False):
	"""Generate scatter plot and gene body snATAC scatter plot using tSNE coordinates from methylation(snmC-seq) data.

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

	genes = genes_query.split()

	gene_name_str = ""
	x, y, text, mch = list(), list(), list(), list()

	if len(genes) == 1:
		points = get_gene_snATAC(ensemble, genes[0], grouping, True, smoothing)
		gene_name = get_gene_by_id([ genes[0] ])[0]['gene_name']
		title = 'Gene body snATAC normalized counts: ' + gene_name
	else:
		points = get_mult_gene_snATAC(ensemble, genes, grouping, smoothing)
		gene_infos = get_gene_by_id(genes)
		for i, gene in enumerate(gene_infos):
			if i > 0 and i % 10 == 0:
				gene_name_str += "<br>"
			gene_name_str += gene['gene_name'] + '+'
		gene_name_str = gene_name_str[:-1]
		title = 'Avg. Gene body snATAC normalized counts: <br>' + gene_name_str

	if points is None:
		raise FailToGraphException

	### TSNE ### 
	if grouping != 'dataset' and grouping != 'target_region':
		if grouping+'_ATAC' not in points.columns: # If no cluster annotations available, group by cluster number instead
			grouping = "cluster"
			if len(genes) == 1:
				points = get_gene_snATAC(ensemble, genes[0], grouping, True, smoothing)
			else:
				points = get_mult_gene_snATAC(ensemble, genes, grouping, smoothing)
			print("**** Grouping by cluster")

	datasets = points['dataset'].unique().tolist()
	annotation_additional_y = 0.00 
	if grouping == 'dataset':
		unique_groups = datasets
		num_clusters = len(unique_groups)
	elif grouping == 'target_region':
		points['target_region'].fillna('N/A', inplace=True)
		unique_groups = points['target_region'].unique().tolist()
		num_clusters = len(unique_groups)
	else:
		if grouping == 'cluster':
			annotation_additional_y = 0.025 # Necessary because legend items overlap with legend title (annotation) when there are many legend items
		num_clusters = points['cluster_ATAC'].max()
		unique_groups = points[grouping+'_ATAC'].unique().tolist()
	
	colors = generate_cluster_colors(len(unique_groups), grouping)
	symbols = ['circle', 'square', 'cross', 'triangle-up', 'triangle-down', 'octagon', 'star', 'diamond']
	
	traces_tsne = OrderedDict()

	legend_x = -.17
	layout_width = 1100;

	grouping_clustering = grouping
	if grouping != 'dataset' and grouping != 'target_region':
		layout_width = 1000;
		legend_x = -.14
		grouping_clustering = grouping+'_ATAC'

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
				   #'symbol': symbols[datasets.index(dataset)],
			},
			hoverinfo='text'))
		trace2d['x'] = points_group['tsne_x_ATAC'].values.tolist()
		trace2d['y'] = points_group['tsne_y_ATAC'].values.tolist()
		# for point in points_group.itertuples(index=False):  # Maybe there's a more elegant way to do this... EAM
		# 	text = OrderedDict([('Cluster', point[4]),('Dataset', point[2]),])
		# 	if point[3]!='Null':
		# 		text['Annotation'] = point[3]
		# 	if point[-1]!='None':
		# 		text['RS2 Target Region'] = point[-1]
		# 	trace2d['text'] = [build_hover_text(OrderedDict(text))]
		trace2d['text'] = [build_hover_text(OrderedDict([('Annotation', point[3]),
														 ('Cluster', point[4]),
														 ('RS2 Target Region', point[-1]),
														 ('Dataset', point[2]),]))
						   for point in points_group.itertuples(index=False)]

	### snATAC normalized counts scatter plot ### 
	x = points['tsne_x_ATAC'].tolist()
	y = points['tsne_y_ATAC'].tolist()
	ATAC_counts = points['normalized_counts'].copy()
	text_ATAC = [build_hover_text(OrderedDict([('Annotation', point[3]),
											   ('Cluster', point[4]),
											   ('RS2 Target Region', point[-1]),
											   ('Dataset', point[2]),
											   ('<b>Normalized Counts</b>', round(point[-2], 5)),]))
				 for point in points.itertuples(index=False)]


	ATAC_dataframe = pd.DataFrame(ATAC_counts)
	start = ATAC_dataframe.dropna().quantile(ptile_start)[0].tolist()
	end = ATAC_dataframe.dropna().quantile(ptile_end).values[0].tolist()
	end = max(end,start+0.01)
	ATAC_colors = [set_color_by_percentile(x, start, end) for x in ATAC_counts]

	colorbar_tickval = list(arange(start, end, (end - start) / 4))
	colorbar_tickval[0] = start
	colorbar_tickval.append(end)
	colorbar_ticktext = [
		str(round(x, num_sigfigs_ticklabels)) for x in arange(start, end, (end - start) / 4)
	]
	colorbar_ticktext[0] = '<' + str(round(start, num_sigfigs_ticklabels))
	colorbar_ticktext.append('>' + str(round(end, num_sigfigs_ticklabels)))

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
		height=550,
		width=layout_width,
		# title=title,
		# titlefont={'color': 'rgba(1,2,2,1)',
		#            'size': 16},
		legend={'x':legend_x,
				'y':0.95,
				'tracegroupgap': 0.5},
		margin={'l': 0,
				'r': 0,
				'b': 30,
				't': 130,},
		xaxis={
			'domain': [0, 0.49],
			'type': 'linear',
			'ticks': '',
			'dtick': 10,
			'tickwidth': 0,
			'showticklabels': False,
			'showline': True,
			'showgrid': False,
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
			'showline': True,
			'showgrid': False,
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
			'showline': True,
			'showgrid': False,
			'side': 'right',
			'zeroline': False,
			'linecolor': 'black',
			'linewidth': 0.5,
			'mirror': False,
			'range':[bottom_y,top_y]
		},
		hovermode='closest',)

	
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
	fig['layout']['annotations'].extend([Annotation(text=grouping.title(),
													x=legend_x+0.05,
													y=1.02 + annotation_additional_y,
													xanchor="left",
													yanchor="top",
													showarrow=False,
													xref="paper",
													yref="paper",
													font={'size': 12,
														  'color': 'gray',})])
	fig['layout']['annotations'].extend([Annotation(text=title,
													x=0.5,
													y=1.3,
													xanchor="center",
													yanchor="top",
													showarrow=False,
													xref="paper",
													yref="paper",
													font={'size': 16,
														  'color': 'black',})])

	return plotly.offline.plot(
		figure_or_data=fig,
		output_type='div',
		show_link=False,
		include_plotlyjs=False)

@cache.memoize(timeout=3600)
def get_snATAC_heatmap(ensemble, grouping, ptile_start, ptile_end, normalize_row, query):
	"""Generate ATAC heatmap comparing multiple genes.

	Arguments:
		ensemble (str): Name of ensemble.
		grouping (str): Variable to group cells by. "cluster" or "annotation"
		ptile_start (float): Lower end of color percentile. [0, 1].
		ptile_end (float): Upper end of color percentile. [0, 1].
		normalize_row (bool): Whether to normalize by each row (gene). 
		query ([str]): Ensembl IDs of genes to display.

	Returns:
		str: HTML generated by Plot.ly
	"""
	
	if normalize_row:
		normal_or_original = '(normalized by gene)'
	else:
		normal_or_original = ''

	title = "Gene body snATAC normalized counts by cluster " + normal_or_original + ":<br>"
	genes = query.split()

	gene_info_df = pd.DataFrame()
	gene_infos = get_gene_by_id(genes)
	for i, gene in enumerate(gene_infos):
		gene_name = gene['gene_name']
		if i > 0 and i % 10 == 0:
			title += "<br>"
		title += gene_name + "+"
		gene_info_df[gene_name] = mean_cluster(get_gene_snATAC(ensemble, gene['gene_id'], grouping, True), grouping, 'ATAC')

	title = title[:-1] # Gets rid of last '+'

	gene_info_df.reset_index(inplace=True)
	if grouping == 'annotation':
		gene_info_df['annotation_cat'] = pd.Categorical(gene_info_df['annotation_ATAC'], cluster_annotation_order)
		gene_info_df.sort_values(by='annotation_cat', inplace=True)
		gene_info_df.drop('annotation_cat', axis=1, inplace=True)
		gene_info_df.set_index(grouping+'_ATAC', inplace=True)
	elif grouping == 'cluster':
		gene_info_df.sort_values(by='cluster_ATAC', inplace=True)
		gene_info_df.set_index(grouping+'_ATAC', inplace=True)
	elif grouping == 'dataset' or grouping == 'target_region':
		gene_info_df.sort_values(by=grouping, inplace=True)
		gene_info_df.set_index(grouping, inplace=True)
	else:
		raise FailToGraphException

	# For some reason, Plotly doesn't allow 'None' as a group on the x-axis for heatmaps.
	if gene_info_df.index.tolist() == ['None']: 
		gene_info_df.index = ['N/A']

	normal_or_original = 'Original'
	if normalize_row:
		for gene in gene_info_df:
			# z-score
			# gene_info_df[gene] = (gene_info_df[gene] - gene_info_df[gene].mean()) / gene_info_df[gene].std()
			# min-max
			gene_range = gene_info_df[gene].max() - gene_info_df[gene].min()
			if (gene_range==0):
				gene_range = 1
			gene_info_df[gene] = (gene_info_df[gene] - gene_info_df[gene].min()) / gene_range
		normal_or_original = 'Normalized'

	gene_info_dict = gene_info_df.to_dict(into=OrderedDict)    

	x, y, text, hover, snATAC_counts = list(), list(), list(), list(), list()
	i = 0
	name_prepend = ""
	if grouping == 'cluster':
		name_prepend = 'cluster_'
	for key in list(gene_info_dict.keys()):
		j = 0
		y.append(key)
		snATAC_counts.append(list(gene_info_dict[key].values()))
		for cluster in list(gene_info_dict[key].keys()):
			x.append(name_prepend+str(cluster))
			text.append(build_hover_text(OrderedDict([('Gene', key),
													  (grouping.title(), x[j]),
													  ('Normalized Counts', snATAC_counts[i][j])
													  ])))
			j += 1
		hover.append(text)
		text = []
		i += 1

	flat_snATAC_counts = list(chain.from_iterable(snATAC_counts))
	snATAC_counts_dataframe = pd.DataFrame(flat_snATAC_counts).dropna()
	start = snATAC_counts_dataframe.quantile(ptile_start).values[0].tolist()
	end = snATAC_counts_dataframe.quantile(ptile_end).values[0].tolist()
	end = max(end,start+0.01)

	colorbar_tickval = list(arange(start, end, (end - start) / 4))
	colorbar_tickval[0] = start
	colorbar_tickval.append(end)
	colorbar_ticktext = [
		str(round(x, num_sigfigs_ticklabels)) for x in arange(start, end, (end - start) / 4)
	]
	if normalize_row == True:
		colorbar_ticktext[0] = str(round(start, num_sigfigs_ticklabels))
	else:
		colorbar_ticktext[0] = '<' + str(round(start, num_sigfigs_ticklabels))
	colorbar_ticktext.append('>' + str(round(end, num_sigfigs_ticklabels)))

	# Due to a weird bug(?) in plotly, the number of elements in tickvals and ticktext 
	# must be greater than or equal to number of genes in query. Else, javascript throws 
	# Uncaught Typeerrors when trying to hover over genes. (Tomo 12/11/17)
	while len(colorbar_tickval) < len(genes):
		colorbar_tickval.insert(0,start)
		if normalize_row == True:
			colorbar_ticktext.insert(0, str(round(start, num_sigfigs_ticklabels)))
		else:
			colorbar_ticktext.insert(0, '<' + str(round(start, num_sigfigs_ticklabels)))

	trace = Heatmap(
		x=x,
		y=y,
		z=snATAC_counts,
		text=hover,
		colorscale='Viridis',
		colorbar={
			'x': 1.0,
			'len': 1,
			'title': 'snATAC normalized counts',
			'titleside': 'right',
			'tickmode': 'array',
			'tickvals': colorbar_tickval,
			'ticktext': colorbar_ticktext,
			'thickness': 10,
			'tickfont': {'size': 10}
			},
		hoverinfo='text',
		zmin=start,zmax=end,zauto=False, # Clip the extreme edges of the colorscale
		)

	layout = Layout(
		autosize=True,
		height=max(600*len(genes)/20,550), # EAM Adjust the height of the heatmap according to the number of genes displayed
		width=1000,
		# title=title,
		# titlefont={'color': 'rgba(1,2,2,1)',
		#            'size': 16},
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
			x=-0.1,
			xanchor='left',
			y=1.43,
			yanchor='top'
		)
	])

	layout['updatemenus'] = updatemenus

	layout['annotations'].extend([Annotation(text=title,
											 x=0.5,
											 y=1.4,
											 xanchor="center",
											 yanchor="top",
											 showarrow=False,
											 xref="paper",
											 yref="paper",
											 font={'size': 16,
												   'color': 'black',})])


	return plotly.offline.plot(
		{
			'data': [trace],
			'layout': layout
		},
		output_type='div',
		show_link=False,
		include_plotlyjs=False)

@cache.memoize(timeout=3600)
def get_snATAC_box(ensemble, gene, grouping, outliers):
	"""Generate gene body mCH box plot.

	Traces are grouped by cluster.

	Arguments:
		ensemble (str): Name of ensemble.
		gene (str):  Ensembl ID of gene for that ensemble.
		grouping (str): Variable to group cells by. "cluster", "annotation".
		outliers (bool): Whether if outliers should be displayed.

	Returns:
		str: HTML generated by Plot.ly.
	"""
	points = get_gene_snATAC(ensemble, gene, grouping, outliers)

	if points is None:
		raise FailToGraphException

	if grouping == "dataset":
		unique_groups = points["dataset"].unique()
	elif grouping == "target_region":
		points['target_region'].fillna('N/A', inplace=True)
		unique_groups = points["target_region"].unique()
	elif grouping == 'annotation' or grouping == 'cluster':
		if grouping == 'annotation' and grouping+'_ATAC' not in points.columns: # If no cluster annotations available, group by cluster number instead
			grouping = "cluster"
			points = get_gene_snATAC(ensemble, gene, grouping, outliers)
			print("**** Grouping by cluster")
		unique_groups = points[grouping+'_ATAC'].unique()
	else: 
		raise FailToGraphException

	num_clusters = len(unique_groups)
	colors = generate_cluster_colors(num_clusters, grouping)

	name_prepend = ""
	x_label = grouping
	if grouping != "dataset" or grouping != "target_region":
		if grouping == "cluster":
			name_prepend="cluster_"
		grouping += "_ATAC"
	if outliers:
		boxpoints='suspectedoutliers';
	else:
		boxpoints=False

	traces = OrderedDict()
	for point in points.to_dict('records'):
		color = colors[int(np.where(unique_groups==point[grouping])[0]) % len(colors)]
		group = point[grouping]
		trace = traces.setdefault(group, Box(
				y=list(),
				name=name_prepend + str(group),
				marker={
					'color': color,
					'outliercolor': color,
					'size': 6
				},
				boxpoints=boxpoints,
				visible=True,
				showlegend=False,
				))
		trace['y'].append(point['normalized_counts'])

	gene_name = get_gene_by_id([ gene ])[0]['gene_name']

	layout = Layout(
		autosize=True,
		height=450,
		width=1000,
		title='Gene body snATAC normalized counts in each cluster:<br>' + gene_name,
		titlefont={'color': 'rgba(1,2,2,1)',
				   'size': 20},
		xaxis={
			'title': x_label.title(),
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
			'title': gene_name + ' snATAC normalized counts',
			'titlefont': {
				'size': 15
			},
			'type': 'linear',
			'anchor': 'x',
			'ticks': 'outside',
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

### RNA
@cache.memoize(timeout=3600)
def get_gene_RNA(ensemble, gene, grouping, outliers):
	"""Return RNA data points for a given gene.

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
	result = db.get_engine(current_app, 'RNA_data').execute("SELECT gene_id FROM genes WHERE gene_id LIKE %s", (gene+"%",)).fetchone()
	gene_table_name = 'gene_' + result['gene_id'].replace('.','_')
	
	query = "SELECT cells.cell_id, cells.cell_name, cells.dataset, \
		%(ensemble)s.annotation_RNA, %(ensemble)s.cluster_RNA, \
		%(ensemble)s.tsne_x_RNA, %(ensemble)s.tsne_y_RNA, \
		%(gene_table_name)s.normalized_counts, \
		datasets.target_region \
		FROM cells \
		INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
		LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
		LEFT JOIN datasets ON cells.dataset = datasets.dataset \
		ORDER BY RAND() LIMIT %(ncells)s" % {'ensemble': ensemble, 
																   'gene_table_name': gene_table_name,
																   'ncells': ncells_max}

	try:
		df = pd.read_sql(query, db.get_engine(current_app, 'RNA_data'))
	except exc.ProgrammingError as e:
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_gene_RNA): {}".format(str(now), e))
		sys.stdout.flush()
		return None

	if df.empty: # If no data in column, return None 
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_gene_RNA): No RNA data for {}".format(str(now), ensemble))
		sys.stdout.flush()
		return None

	if grouping == 'annotation':
		df.fillna({'annotation_RNA': 'None'}, inplace=True)
		df['annotation_cat'] = pd.Categorical(df['annotation_RNA'], cluster_annotation_order)
		df.sort_values(by='annotation_cat', inplace=True)
		df.drop('annotation_cat', axis=1, inplace=True)
	elif grouping == 'cluster':
		df.sort_values(by='cluster_RNA', inplace=True)

	df['normalized_counts'].fillna(0, inplace=True)
	
	return df

@cache.memoize(timeout=1800)
def get_mult_gene_RNA(ensemble, genes, grouping):
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
			%(ensemble)s.annotation_RNA, %(ensemble)s.cluster_RNA, \
			%(ensemble)s.tsne_x_RNA, %(ensemble)s.tsne_y_RNA, \
			%(gene_table_name)s.normalized_counts, \
			datasets.target_region \
			FROM cells \
			INNER JOIN %(ensemble)s ON cells.cell_id = %(ensemble)s.cell_id \
			LEFT JOIN %(gene_table_name)s ON %(ensemble)s.cell_id = %(gene_table_name)s.cell_id \
			LEFT JOIN datasets ON cells.dataset = datasets.dataset \
			ORDER BY RAND() LIMIT %(ncells)s" % {'ensemble': ensemble, 
																	   'gene_table_name': gene_table_name,
																	   'ncells': ncells_max}
		try:
			df_all = df_all.append(pd.read_sql(query, db.get_engine(current_app, 'RNA_data')))
		except exc.ProgrammingError as e:
			now = datetime.datetime.now()
			print("[{}] ERROR in app(get_mult_gene_RNA): {}".format(str(now), e))
			sys.stdout.flush()
			return None
		
		if first:
			df_coords = df_all
		first = False

	if df_all.empty: # If no data in column, return None 
		now = datetime.datetime.now()
		print("[{}] ERROR in app(get_gene_RNA): No RNA data for {}".format(str(now), ensemble))
		sys.stdout.flush()
		return None

	df_all['normalized_counts'].fillna(0, inplace=True)

	df_avg_methylation = df_all.groupby(by='cell_id', as_index=False)['normalized_counts'].mean()
	df_coords.update(df_avg_methylation)

	if grouping == 'annotation':
		df_coords.fillna({'annotation_RNA': 'None'}, inplace=True)
		df_coords['annotation_cat'] = pd.Categorical(df_coords['annotation_RNA'], cluster_annotation_order)
		df_coords.sort_values(by='annotation_cat', inplace=True)
		df_coords.drop('annotation_cat', axis=1, inplace=True)
	elif grouping == 'cluster':
		df_coords.sort_values(by='cluster_RNA', inplace=True)
	return df_coords

@cache.memoize(timeout=1800)
def get_RNA_scatter(ensemble, genes_query, grouping, ptile_start, ptile_end, tsne_outlier_bool):
	"""Generate RNA scatter plot using tSNE coordinates from methylation(snmC-seq) data.

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

	genes = genes_query.split()

	gene_name_str = ""
	x, y, text, mch = list(), list(), list(), list()

	if len(genes) == 1:
		points = get_gene_RNA(ensemble, genes[0], grouping, True)
		gene_name = get_gene_by_id([ genes[0] ])[0]['gene_name']
		title = 'Gene body RNA normalized counts: ' + gene_name
	else:
		points = get_mult_gene_RNA(ensemble, genes, grouping)
		gene_infos = get_gene_by_id(genes)
		for i, gene in enumerate(gene_infos):
			if i > 0 and i % 10 == 0:
				gene_name_str += "<br>"
			gene_name_str += gene['gene_name'] + '+'
		gene_name_str = gene_name_str[:-1]
		title = 'Avg. Gene body RNA normalized counts: <br>' + gene_name_str

	if points is None:
		raise FailToGraphException

	### TSNE ### 
	if grouping != 'dataset' and grouping != 'target_region':
		if grouping+'_RNA' not in points.columns: # If no cluster annotations available, group by cluster number instead
			grouping = "cluster"
			if len(genes) == 1:
				points = get_gene_RNA(ensemble, genes[0], grouping, True)
			else:
				points = get_mult_gene_RNA(ensemble, genes, grouping)
			print("**** Grouping by cluster")

	datasets = points['dataset'].unique().tolist()
	annotation_additional_y = 0.00 
	if grouping == 'dataset':
		unique_groups = datasets
		num_clusters = len(unique_groups)
	elif grouping == 'target_region':
		points['target_region'].fillna('N/A', inplace=True)
		unique_groups = points['target_region'].unique().tolist()
		num_clusters = len(unique_groups)
	else:
		if grouping == 'cluster':
			annotation_additional_y = 0.025 # Necessary because legend items overlap with legend title (annotation) when there are many legend items
		num_clusters = points['cluster_RNA'].max()
		unique_groups = points[grouping+'_RNA'].unique().tolist()
	
	colors = generate_cluster_colors(len(unique_groups), grouping)
	symbols = ['circle', 'square', 'cross', 'triangle-up', 'triangle-down', 'octagon', 'star', 'diamond']
	
	traces_tsne = OrderedDict()

	legend_x = -.17
	layout_width = 1100;

	grouping_clustering = grouping
	if grouping != 'dataset' and grouping != 'target_region':
		layout_width = 1000;
		legend_x = -.14
		grouping_clustering = grouping+'_RNA'

	if tsne_outlier_bool:
		top_x = points['tsne_x_RNA'].quantile(0.999)
		bottom_x = points['tsne_x_RNA'].quantile(0.001) 
		top_y = points['tsne_y_RNA'].quantile(0.999)
		bottom_y = points['tsne_y_RNA'].quantile(0.001) 
		
	else:
		top_x = points['tsne_x_RNA'].max()
		bottom_x = points['tsne_x_RNA'].min()
		top_y = points['tsne_y_RNA'].max()
		bottom_y = points['tsne_y_RNA'].min()

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
				   #'symbol': symbols[datasets.index(dataset)],
			},
			hoverinfo='text'))
		trace2d['x'] = points_group['tsne_x_RNA'].values.tolist()
		trace2d['y'] = points_group['tsne_y_RNA'].values.tolist()
		trace2d['text'] = [build_hover_text(OrderedDict([('Annotation', point[3]),
														 ('Cluster', point[4]),
														 ('RS2 Target Region', point[-1]),
														 ('Dataset', point[2]),]))
						   for point in points_group.itertuples(index=False)]

	### RNA normalized counts scatter plot ### 
	x = points['tsne_x_RNA'].tolist()
	y = points['tsne_y_RNA'].tolist()
	RNA_counts = points['normalized_counts'].copy()
	text_RNA = [build_hover_text(OrderedDict([('Annotation', point[3]),
											   ('Cluster', point[4]),
											   ('RS2 Target Region', point[-1]),
											   ('Dataset', point[2]),
											   ('<b>Normalized Counts</b>', round(point[-2], 5)),]))
				 for point in points.itertuples(index=False)]


	RNA_dataframe = pd.DataFrame(RNA_counts)
	start = RNA_dataframe.dropna().quantile(ptile_start)[0].tolist()
	end = RNA_dataframe.dropna().quantile(ptile_end).values[0].tolist()
	end = max(end,start+0.01)
	RNA_colors = [set_color_by_percentile(x, start, end) for x in RNA_counts]

	colorbar_tickval = list(arange(start, end, (end - start) / 4))
	colorbar_tickval[0] = start
	colorbar_tickval.append(end)
	colorbar_ticktext = [
		str(round(x, num_sigfigs_ticklabels)) for x in arange(start, end, (end - start) / 4)
	]
	colorbar_ticktext[0] = '<' + str(round(start, num_sigfigs_ticklabels))
	colorbar_ticktext.append('>' + str(round(end, num_sigfigs_ticklabels)))

	trace_RNA = Scatter(
		mode='markers',
		x=x,
		y=y,
		text=text_RNA,
		marker={
			'color': RNA_colors,
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
		height=550,
		width=layout_width,
		# title=title,
		# titlefont={'color': 'rgba(1,2,2,1)',
		#            'size': 16},
		legend={'x':legend_x,
				'y':0.95,
				'tracegroupgap': 0.5},
		margin={'l': 0,
				'r': 0,
				'b': 30,
				't': 130,},
		xaxis={
			'domain': [0, 0.49],
			'type': 'linear',
			'ticks': '',
			'dtick': 10,
			'tickwidth': 0,
			'showticklabels': False,
			'showline': True,
			'showgrid': False,
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
			'showline': True,
			'showgrid': False,
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
			'showline': True,
			'showgrid': False,
			'side': 'right',
			'zeroline': False,
			'linecolor': 'black',
			'linewidth': 0.5,
			'mirror': False,
			'range':[bottom_y,top_y]
		},
		hovermode='closest',)

	
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
	fig.append_trace(trace_RNA, 1,2)

	fig['layout'].update(layout)
	fig['layout']['annotations'].extend([Annotation(text=grouping.title(),
													x=legend_x+0.05,
													y=1.02 + annotation_additional_y,
													xanchor="left",
													yanchor="top",
													showarrow=False,
													xref="paper",
													yref="paper",
													font={'size': 12,
														  'color': 'gray',})])
	fig['layout']['annotations'].extend([Annotation(text=title,
													x=0.5,
													y=1.3,
													xanchor="center",
													yanchor="top",
													showarrow=False,
													xref="paper",
													yref="paper",
													font={'size': 16,
														  'color': 'black',})])

	return plotly.offline.plot(
		figure_or_data=fig,
		output_type='div',
		show_link=False,
		include_plotlyjs=False)

@cache.memoize(timeout=3600)
def get_RNA_heatmap(ensemble, grouping, ptile_start, ptile_end, normalize_row, query):
	"""Generate RNA heatmap comparing multiple genes.

	Arguments:
		ensemble (str): Name of ensemble.
		grouping (str): Variable to group cells by. "cluster" or "annotation"
		ptile_start (float): Lower end of color percentile. [0, 1].
		ptile_end (float): Upper end of color percentile. [0, 1].
		normalize_row (bool): Whether to normalize by each row (gene). 
		query ([str]): Ensembl IDs of genes to display.

	Returns:
		str: HTML generated by Plot.ly
	"""
	
	if normalize_row:
		normal_or_original = '(normalized by gene)'
	else:
		normal_or_original = ''

	title = "Gene body RNA normalized counts by cluster " + normal_or_original + ":<br>"
	genes = query.split()

	gene_info_df = pd.DataFrame()
	gene_infos = get_gene_by_id(genes)
	for i, gene in enumerate(gene_infos):
		gene_name = gene['gene_name']
		if i > 0 and i % 10 == 0:
			title += "<br>"
		title += gene_name + "+"
		gene_info_df[gene_name] = mean_cluster(get_gene_RNA(ensemble, gene['gene_id'], grouping, True), grouping, 'RNA')

	title = title[:-1] # Gets rid of last '+'

	gene_info_df.reset_index(inplace=True)
	if grouping == 'annotation':
		gene_info_df['annotation_cat'] = pd.Categorical(gene_info_df['annotation_RNA'], cluster_annotation_order)
		gene_info_df.sort_values(by='annotation_cat', inplace=True)
		gene_info_df.drop('annotation_cat', axis=1, inplace=True)
		gene_info_df.set_index(grouping+'_RNA', inplace=True)
	elif grouping == 'cluster':
		gene_info_df.sort_values(by='cluster_RNA', inplace=True)
		gene_info_df.set_index(grouping+'_RNA', inplace=True)
	elif grouping == 'dataset' or grouping == 'target_region':
		gene_info_df.sort_values(by=grouping, inplace=True)
		gene_info_df.set_index(grouping, inplace=True)
	else:
		raise FailToGraphException

	# For some reason, Plotly doesn't allow 'None' as a group on the x-axis for heatmaps.
	if gene_info_df.index.tolist() == ['None']: 
		gene_info_df.index = ['N/A']

	normal_or_original = 'Original'
	if normalize_row:
		for gene in gene_info_df:
			# z-score
			# gene_info_df[gene] = (gene_info_df[gene] - gene_info_df[gene].mean()) / gene_info_df[gene].std()
			# min-max
			gene_range = gene_info_df[gene].max() - gene_info_df[gene].min()
			if (gene_range==0):
				gene_range = 1
			gene_info_df[gene] = (gene_info_df[gene] - gene_info_df[gene].min()) / gene_range
		normal_or_original = 'Normalized'

	gene_info_dict = gene_info_df.to_dict(into=OrderedDict)    

	x, y, text, hover, RNA_counts = list(), list(), list(), list(), list()
	i = 0
	name_prepend = ""
	if grouping == 'cluster':
		name_prepend = 'cluster_'
	for key in list(gene_info_dict.keys()):
		j = 0
		y.append(key)
		RNA_counts.append(list(gene_info_dict[key].values()))
		for cluster in list(gene_info_dict[key].keys()):
			x.append(name_prepend+str(cluster))
			text.append(build_hover_text(OrderedDict([('Gene', key),
													  (grouping.title(), x[j]),
													  ('Normalized Counts', RNA_counts[i][j])
													  ])))
			j += 1
		hover.append(text)
		text = []
		i += 1

	flat_RNA_counts = list(chain.from_iterable(RNA_counts))
	RNA_counts_dataframe = pd.DataFrame(flat_RNA_counts).dropna()
	start = RNA_counts_dataframe.quantile(ptile_start).values[0].tolist()
	end = RNA_counts_dataframe.quantile(ptile_end).values[0].tolist()
	end = max(end,start+0.01)

	colorbar_tickval = list(arange(start, end, (end - start) / 4))
	colorbar_tickval[0] = start
	colorbar_tickval.append(end)
	colorbar_ticktext = [
		str(round(x, num_sigfigs_ticklabels)) for x in arange(start, end, (end - start) / 4)
	]
	if normalize_row == True:
		colorbar_ticktext[0] = str(round(start, num_sigfigs_ticklabels))
	else:
		colorbar_ticktext[0] = '<' + str(round(start, num_sigfigs_ticklabels))
	colorbar_ticktext.append('>' + str(round(end, num_sigfigs_ticklabels)))

	# Due to a weird bug(?) in plotly, the number of elements in tickvals and ticktext 
	# must be greater than or equal to number of genes in query. Else, javascript throws 
	# Uncaught Typeerrors when trying to hover over genes. (Tomo 12/11/17)
	while len(colorbar_tickval) < len(genes):
		colorbar_tickval.insert(0,start)
		if normalize_row == True:
			colorbar_ticktext.insert(0, str(round(start, num_sigfigs_ticklabels)))
		else:
			colorbar_ticktext.insert(0, '<' + str(round(start, num_sigfigs_ticklabels)))

	trace = Heatmap(
		x=x,
		y=y,
		z=RNA_counts,
		text=hover,
		colorscale='Viridis',
		colorbar={
			'x': 1.0,
			'len': 1,
			'title': 'RNA normalized counts',
			'titleside': 'right',
			'tickmode': 'array',
			'tickvals': colorbar_tickval,
			'ticktext': colorbar_ticktext,
			'thickness': 10,
			'tickfont': {'size': 10}
			},
		hoverinfo='text',
		zmin=start,zmax=end,zauto=False, # Clip the extreme edges of the colorscale
		)

	layout = Layout(
		autosize=True,
		height=max(600*len(genes)/20,550), # EAM Adjust the height of the heatmap according to the number of genes displayed
		width=1000,
		# title=title,
		# titlefont={'color': 'rgba(1,2,2,1)',
		#            'size': 16},
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
			x=-0.1,
			xanchor='left',
			y=1.43,
			yanchor='top'
		)
	])

	layout['updatemenus'] = updatemenus

	layout['annotations'].extend([Annotation(text=title,
											 x=0.5,
											 y=1.4,
											 xanchor="center",
											 yanchor="top",
											 showarrow=False,
											 xref="paper",
											 yref="paper",
											 font={'size': 16,
												   'color': 'black',})])


	return plotly.offline.plot(
		{
			'data': [trace],
			'layout': layout
		},
		output_type='div',
		show_link=False,
		include_plotlyjs=False)

@cache.memoize(timeout=3600)
def get_RNA_box(ensemble, gene, grouping, outliers):
	"""Generate gene body mCH box plot.

	Traces are grouped by cluster.

	Arguments:
		ensemble (str): Name of ensemble.
		gene (str):  Ensembl ID of gene for that ensemble.
		grouping (str): Variable to group cells by. "cluster", "annotation".
		outliers (bool): Whether if outliers should be displayed.

	Returns:
		str: HTML generated by Plot.ly.
	"""
	points = get_gene_RNA(ensemble, gene, grouping, outliers)

	if points is None:
		raise FailToGraphException

	if grouping == "dataset":
		unique_groups = points["dataset"].unique()
	elif grouping == "target_region":
		points['target_region'].fillna('N/A', inplace=True)
		unique_groups = points["target_region"].unique()
	elif grouping == 'annotation' or grouping == 'cluster':
		if grouping == 'annotation' and grouping+'_RNA' not in points.columns: # If no cluster annotations available, group by cluster number instead
			grouping = "cluster"
			points = get_gene_RNA(ensemble, gene, grouping, outliers)
			print("**** Grouping by cluster")
		unique_groups = points[grouping+'_RNA'].unique()
	else: 
		raise FailToGraphException

	num_clusters = len(unique_groups)
	colors = generate_cluster_colors(num_clusters, grouping)

	name_prepend = ""
	x_label = grouping
	if grouping != "dataset" or grouping != "target_region":
		if grouping == "cluster":
			name_prepend="cluster_"
		grouping += "_RNA"
	if outliers:
		boxpoints='suspectedoutliers';
	else:
		boxpoints=False

	traces = OrderedDict()
	for point in points.to_dict('records'):
		color = colors[int(np.where(unique_groups==point[grouping])[0]) % len(colors)]
		group = point[grouping]
		trace = traces.setdefault(group, Box(
				y=list(),
				name=name_prepend + str(group),
				marker={
					'color': color,
					'outliercolor': color,
					'size': 6
				},
				boxpoints=boxpoints,
				visible=True,
				showlegend=False,
				))
		trace['y'].append(point['normalized_counts'])

	gene_name = get_gene_by_id([ gene ])[0]['gene_name']

	layout = Layout(
		autosize=True,
		height=450,
		width=1000,
		title='Gene body RNA normalized counts in each cluster:<br>' + gene_name,
		titlefont={'color': 'rgba(1,2,2,1)',
				   'size': 20},
		xaxis={
			'title': x_label.title(),
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
			'title': gene_name + ' RNA normalized counts',
			'titlefont': {
				'size': 15
			},
			'type': 'linear',
			'anchor': 'x',
			'ticks': 'outside',
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

