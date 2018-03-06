# App-specific configuration

# Serve the Bootstrap template locally
# Set to False due to SCCF's subfolder reverse proxy shenanigans, 
# which will cause HTML paths to fail
BOOTSTRAP_SERVE_LOCAL = False

# Disable pretty-printing to conserve network transfer
JSONIFY_PRETTYPRINT_REGULAR = False

# Minify HTML to conserve network transfer
MINIFY_PAGE = True  

MAIL_SERVER = 'smtp.sendgrid.net'
MAIL_PORT = 465
MAIL_USE_TLS = False
MAIL_USE_SSL = True
MAIL_DEBUG =  True
MAIL_USERNAME = 'brainome-admin@ucsd.edu'
MAIL_PASSWORD = ''

MYSQL_USER = ''
MYSQL_PW = ''
MYSQL_DB = ''
SQLALCHEMY_BINDS = { 'data': 'mysql://' + MYSQL_USER + ':' + MYSQL_PW + '@banjo/' + MYSQL_DB }

DATA_DIR = ''
PUBLISHED_DATA_DIR = ''
ALL_DATA_DIR = ''

SECRET_KEY = ''
