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
from .assets import app_css, app_js, vendor_css, vendor_js
import urllib.parse
from flask_wtf import CsrfProtect


# Necessary because brainome doesn't have mysql installed
import pymysql
pymysql.install_as_MySQLdb()

cache = Cache(config={'CACHE_TYPE': 'simple', 'CACHE_THRESHOLD': 1000})
nav = Nav()
mail = Mail()
db = SQLAlchemy()
csrf = CsrfProtect()
compress = Compress()
htmlmin = HTMLMIN()
basedir = os.path.abspath(os.path.dirname(__file__))

# Set up Flask-Login
login_manager = LoginManager()
login_manager.session_protection = 'strong'
login_manager.login_view = 'frontend.login'

class MiniJSONEncoder(JSONEncoder):
    """Minify JSON output."""
    item_separator = ','
    key_separator = ':'

def create_app(configfile=None):
    app = Flask(__name__)
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

    app.config['MAIL_SERVER'] = 'smtp.sendgrid.net'
    app.config['MAIL_PORT'] = 465
    app.config['MAIL_USE_TLS'] = False
    app.config['MAIL_USE_SSL'] = True
    app.config['MAIL_DEBUG'] = True
    app.config['MAIL_USERNAME'] = ''
    app.config['MAIL_PASSWORD'] = ''
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

    with app.app_context():
        from .frontend import frontend
        app.register_blueprint(frontend)

        from .content import content
        app.register_blueprint(content)

    app.json_encoder = MiniJSONEncoder

    nav.init_app(app)
    mail.init_app(app)
    csrf.init_app(app)
    db.init_app(app)
    login_manager.init_app(app)
    compress.init_app(app)
    htmlmin.init_app(app)
    RQ(app)

    return app

if __name__ == '__main__':
  # app.run(debug=True)
  create_app().start(debug=True)

