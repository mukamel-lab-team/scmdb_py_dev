import os
from flask import Flask, current_app
from flask_mail import Mail
from flask_appconfig import AppConfig
from flask_bootstrap import Bootstrap
from flask.ext.cache import Cache
from flask_nav import Nav
from flask_assets import Environment
from flask_compress import Compress
from flask_htmlmin import HTMLMIN
from flask.json import JSONEncoder
from flask_login import LoginManager
from flask_sqlalchemy import SQLAlchemy
from flask_rq import RQ
from flask_jsglue import JSGlue
from .assets import app_css, app_js, vendor_css, vendor_js, browser_js, browser_css, tabular_rs1_js, tabular_rs2_js, tabular_ensemble_js, tabular_css, request_new_ensemble_js
import urllib.parse
from flask_wtf import CsrfProtect


# Necessary because brainome doesn't have mysql installed
import pymysql
pymysql.install_as_MySQLdb()

class MiniJSONEncoder(JSONEncoder):
    """Minify JSON output."""
    item_separator = ','
    key_separator = ':'

cache = Cache(config={'CACHE_TYPE': 'simple', 'CACHE_THRESHOLD': 1000})
nav = Nav()
mail = Mail()
db = SQLAlchemy()
csrf = CsrfProtect()
compress = Compress()
htmlmin = HTMLMIN()
jsglue = JSGlue()
basedir = os.path.abspath(os.path.dirname(__file__))

# Set up Flask-Login
login_manager = LoginManager()
login_manager.session_protection = 'strong'
login_manager.login_view = 'frontend.login'

#app = Flask(__name__, template_folder=os.path.join(os.path.dirname(__file__), '/templates'))
app = Flask(__name__)


#def create_app(configfile=None):
AppConfig(app)
Bootstrap(app)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///' + os.path.join(basedir, 'user-login.sqlite')
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
REDIS_URL = 'http://localhost:6379'
urllib.parse.uses_netloc.append('redis')
url = urllib.parse.urlparse(REDIS_URL)
app.config['RQ_DEFAULT_HOST'] = url.hostname
app.config['RQ_DEFAULT_PORT'] = url.port
app.config['RQ_DEFAULT_PASSWORD'] = url.password
app.config['RQ_DEFAULT_DB'] = 0

# EAM : Set limit on the number of items in cache (RAM)
cache.init_app(app)

# Set up asset pipeline
assets_env = Environment(app)
dirs = ['assets/styles', 'assets/scripts']
for path in dirs:
    assets_env.append_path(os.path.join(basedir, path))
assets_env.url_expire = True

assets_env.register('app_css', app_css)
assets_env.register('app_js', app_js)
assets_env.register('vendor_css', vendor_css)
assets_env.register('vendor_js', vendor_js)
assets_env.register('browser_js', browser_js)
assets_env.register('browser_css', browser_css)
assets_env.register('tabular_ensemble_js', tabular_ensemble_js)
assets_env.register('tabular_rs1_js', tabular_rs1_js)
assets_env.register('tabular_rs2_js', tabular_rs2_js)
assets_env.register('request_new_ensemble_js', request_new_ensemble_js)
assets_env.register('tabular_css', tabular_css)

with app.app_context():
    from .frontend import frontend
    #app.register_blueprint(frontend, url_prefix="/portal")
    app.register_blueprint(frontend)

    from .content import content
    app.register_blueprint(content)
    #app.register_blueprint(content, url_prefix="/portal")

app.json_encoder = MiniJSONEncoder

nav.init_app(app)
mail.init_app(app)
csrf.init_app(app)
db.init_app(app)
login_manager.init_app(app)
compress.init_app(app)
htmlmin.init_app(app)
jsglue.init_app(app)
RQ(app)


if __name__ == '__main__':
  # app.run(debug=True)
  app.run()

