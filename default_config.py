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

MYSQL_USER = ''
MYSQL_PW = ''
MYSQL_DB_methylation = ''
MYSQL_DB_snATAC = ''
SQLALCHEMY_BINDS = {'methylation_data': 'mysql://' + MYSQL_USER + ':' + MYSQL_PW + '@localhost/' + MYSQL_DB_methylation,
                    'snATAC_data': 'mysql://' + MYSQL_USER + ':' + MYSQL_PW + '@localhost/' + MYSQL_DB_snATAC}

DATA_DIR = ''
PUBLISHED_DATA_DIR = ''
ALL_DATA_DIR = ''

SECRET_KEY = ''
