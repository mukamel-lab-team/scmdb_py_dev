# App-specific configuration

# Serve the Bootstrap template locally.
# This is set to False due to SCCF's subfolder reverse proxy shenanigans, which will cause HTML paths to fail.
BOOTSTRAP_SERVE_LOCAL = False

JSONIFY_PRETTYPRINT_REGULAR = False  # Disable pretty-printing to conserve network transfer.
MINIFY_PAGE = True  # Minify HTML to conserve network transfer.

DATA_DIR = '/srv/scmdb_py_newdata/data/'
