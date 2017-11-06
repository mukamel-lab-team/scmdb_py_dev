"""Defines the web routes of the website.

For actual content generation see the content.py module.
"""
from flask import Blueprint, render_template, jsonify, request, redirect, current_app
from flask_nav.elements import Navbar, Link

from .content import get_cluster_plot, search_gene_names, \
    get_methylation_scatter, get_mch_box, get_mch_box_two_species, \
    find_orthologs, FailToGraphException, get_corr_genes, \
    gene_id_to_name, randomize_cluster_colors, get_mch_heatmap, all_gene_modules, \
    get_genes_of_module
from .nav import nav
from .cache import cache
from os import walk

frontend = Blueprint('frontend', __name__)  # Flask "bootstrap"

# Find all the samples in the data directory
dir_list = next(walk(current_app.config['DATA_DIR']))[1]

dir_list_links = [Link(x, x) for x in dir_list]

nav.register_element('frontend_top',
                     Navbar('', *dir_list_links))


# Visitor routes
@frontend.route('/')
def index():
    # Index is not needed since this site is embedded as a frame.
    # We use a JavaScript redirection here, since a reverse proxy will be confused about subdirectories.
    return 'To be redirected manually, click <a href="./human_combined">here</a>.' + \
           '<script>window.location = "./human_combined"; window.location.replace("./human_combined");</script>'


@frontend.route('/<species>')
def species(species):
    return render_template('speciesview.html', species=species)


@frontend.route('/standalone/<species>/<gene>')
def standalone(species, gene):  # View gene body mCH plots alone.
    return render_template('mch_standalone.html', species=species, gene=gene)


@frontend.route('/compare/<mmu_gid>/<hsa_gid>')
def compare(mmu_gid, hsa_gid):
    return render_template('compareview.html', mmu_gid=mmu_gid, hsa_gid=hsa_gid)


@frontend.route('/box_combined/<mmu_gid>/<hsa_gid>')
def box_combined(mmu_gid, hsa_gid):
    return render_template(
        'combined_box_standalone.html', mmu_gid=mmu_gid, hsa_gid=hsa_gid)


# API routes
@cache.cached(timeout=3600)
@frontend.route('/plot/cluster/<species>/<grouping>')
def plot_cluster(species, grouping):
    try:
        return jsonify(get_cluster_plot(species, grouping))
    except FailToGraphException:
        return 'Failed to produce cluster plot. Contact maintainer.'


@cache.cached(timeout=3600)
@frontend.route('/plot/scatter/<species>/<methylationType>/<level>/<ptile_start>/<ptile_end>')
def plot_methylation_scatter(species, methylationType, level, ptile_start, ptile_end):
    genes = request.args.get('q', 'MustHaveAQueryString')
    try:
        return get_methylation_scatter(species,
                                       methylationType,
                                       genes, level,
                                       float(ptile_start), float(ptile_end))
    except (FailToGraphException, ValueError) as e:
        print(e)
        return 'Failed to produce mCH levels scatter plot. Contact maintainer.'


@cache.cached(timeout=3600)
@frontend.route('/plot/box/<species>/<methylationType>/<gene>/<level>/<outliers_toggle>')
def plot_mch_box(species, methylationType, gene, level, outliers_toggle):
    if outliers_toggle == 'outliers':
        outliers = True
    else:
        outliers = False

    try:
        return get_mch_box(species, methylationType, gene, level, outliers)
    except (FailToGraphException, ValueError) as e:
        print(e)
        return 'Failed to produce mCH levels box plot. Contact maintainer.'


@cache.cached(timeout=3600)
@frontend.route('/plot/box_combined/<species>/<methylationType>/<gene_mmu>/<gene_hsa>/<level>/<outliers_toggle>')
def plot_mch_box_two_species(species, methylationType, gene_mmu, gene_hsa, level, outliers_toggle):
    if outliers_toggle == 'outliers':
        outliers = True
    else:
        outliers = False

    try:
        return get_mch_box_two_species(species, methylationType, gene_mmu, gene_hsa, level, outliers)
    except (FailToGraphException, ValueError) as e:
        print(e)
        return 'Failed to produce mCH levels box plot. Contact maintainer.'


@cache.cached(timeout=3600)
@frontend.route('/gene/names/<species>')
def search_gene_by_name(species):
    query = request.args.get('q', 'MustHaveAQueryString')
    if query == 'none' or query == '':
        return jsonify([])
    else:
        return jsonify(search_gene_names(species, query))


@cache.cached(timeout=3600)
@frontend.route('/gene/id/<species>')
def search_gene_by_id(species):
    query = request.args.get('q', 'MustHaveAQueryString')
    if query == 'none' or query == '':
        return jsonify([])
    else:
        return jsonify(gene_id_to_name(species, query))


@cache.cached(timeout=3600)
@frontend.route('/gene/modules')
def gene_modules():
    query = request.args.get('q')
    if query == None or query == '':
        return jsonify(all_gene_modules())
    else:
        return jsonify(get_genes_of_module(query))

@frontend.route('/gene/orthologs/<species>/<geneID>')
def orthologs(species, geneID):
    geneID = geneID.split('.')[0]
    if species == 'mmu':
        return jsonify(find_orthologs(mmu_gid=geneID))
    else:
        return jsonify(find_orthologs(hsa_gid=geneID))


@frontend.route('/gene/corr/<species>/<geneID>')
def correlated_genes(species, geneID):
    return jsonify(get_corr_genes(species, geneID))


@frontend.route('/plot/randomize_colors')
def randomize_colors():
    return jsonify(randomize_cluster_colors())


@frontend.route('/plot/heat/<species>/<methylationType>/<level>/<ptile_start>/<ptile_end>')
def plot_mch_heatmap(species, methylationType, level, ptile_start, ptile_end):
    query = request.args.get('q', 'MustHaveAQueryString')
    return get_mch_heatmap(species, methylationType, level, ptile_start, ptile_end, query)


@frontend.route('/help')
def help_page():
    return render_template('help.html')
