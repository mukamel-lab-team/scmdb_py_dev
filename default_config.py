# App-specific configuration

# Serve the Bootstrap template locally
# Set to False due to SCCF's subfolder reverse proxy shenanigans, 
# which will cause HTML paths to fail
BOOTSTRAP_SERVE_LOCAL = False

# Disable pretty-printing to conserve network transfer
JSONIFY_PRETTYPRINT_REGULAR = False

# Minify HTML to conserve network transfer
MINIFY_PAGE = True  

# Note: This now points to brainome_210 which has a completely different layout
# DATA_DIR = '/srv/scmdb_py/data'
# PUBLISHED_DATA_DIR = '/srv/scmdb_py/data'
# ALL_DATA_DIR = '/srv/scmdb_py_newdata/data'

DATA_DIR = '/home/cogsci-cnddcollab/brainome_210/data/published_data'
PUBLISHED_DATA_DIR = '/home/cogsci-cnddcollab/brainome_210/data/published_data'
ALL_DATA_DIR = '/home/cogsci-cnddcollab/brainome_210/data/all_data'
