# App-specific configuration

# Serve the Bootstrap template locally
# Set to False due to SCCF's subfolder reverse proxy shenanigans, 
# which will cause HTML paths to fail
BOOTSTRAP_SERVE_LOCAL = False

# Disable pretty-printing to conserve network transfer
JSONIFY_PRETTYPRINT_REGULAR = False

# Minify HTML to conserve network transfer
MINIFY_PAGE = True  

# Statement for enabling the development environment
#DEBUG = True

THREADS_PER_PAGE = 2

MAIL_SERVER = 'smtp.ucsd.edu'
MAIL_PORT = 465
MAIL_USE_TLS = False
MAIL_USE_SSL = True
MAIL_DEBUG =  True
MAIL_USERNAME = 'lab@brainome.ucsd.edu'

REQUEST_EMAIL = 'f7xie@ad.ucsd.edu'

MYSQL_USER = 'cndd_annoj'
MYSQL_PW = 'jonna_ddnc'
MYSQL_DB_methylation = 'CEMBA'
MYSQL_DB_snATAC = 'CEMBA_snATAC'
SQLALCHEMY_BINDS = {'methylation_data': 'mysql://' + MYSQL_USER + ':' + MYSQL_PW + '@banjo/' + MYSQL_DB_methylation,
                    'snATAC_data': 'mysql://' + MYSQL_USER + ':' + MYSQL_PW + '@banjo/' + MYSQL_DB_snATAC}

DATA_DIR = ''
PUBLISHED_DATA_DIR = ''
ALL_DATA_DIR = ''

# Enable protection agains *Cross-site Request Forgery (CSRF)*
CSRF_ENABLED = True

# Use a secure, unique and absolutely secret key for
# signing the data. 
CSRF_SESSION_KEY = "s3cr3t-brainome-csrf"

SECRET_KEY = 's3cr3t'

APPLICATION_ROOT = '/portal'
