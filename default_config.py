# App-specific configuration

# Serve the Bootstrap template locally
# Set to False due to SCCF's subfolder reverse proxy shenanigans, 
# which will cause HTML paths to fail
BOOTSTRAP_SERVE_LOCAL = False

# Disable pretty-printing to conserve network transfer
JSONIFY_PRETTYPRINT_REGULAR = False

# Minify HTML to conserve network transfer
MINIFY_PAGE = True  

MAIL_SERVER = 'smtp.ucsd.edu'
MAIL_PORT = 465
MAIL_USE_TLS = False
MAIL_USE_SSL = True
MAIL_DEBUG =  True
MAIL_USERNAME = 'lab@brainome.ucsd.edu'

REQUEST_EMAIL = 'tishihar94@gmail.com'

MYSQL_USER = 'tishihar'
MYSQL_PW = 'password1234'
MYSQL_DB_methylation = 'CEMBA'
MYSQL_DB_snATAC = 'CEMBA_snATAC'
SQLALCHEMY_BINDS = {'methylation_data': 'mysql://' + MYSQL_USER + ':' + MYSQL_PW + '@localhost/' + MYSQL_DB_methylation,
                    'snATAC_data': 'mysql://' + MYSQL_USER + ':' + MYSQL_PW + '@localhost/' + MYSQL_DB_snATAC}

DATA_DIR = ''
PUBLISHED_DATA_DIR = ''
ALL_DATA_DIR = ''

SECRET_KEY = 's3cr3t'
