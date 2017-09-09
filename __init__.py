# -*- coding: utf-8 -*-
# Main entry point.
#
# In its parent directory, execute:
# $ flask --app=scmdb_py dev

from flask import Flask
from flask_appconfig import AppConfig
from flask_bootstrap import Bootstrap

from .frontend import frontend
from .nav import nav
from .cache import cache
from .compress import compress, htmlmin
from .json import MiniJSONEncoder


def create_app(configfile=None):
    app = Flask(__name__)

    AppConfig(app)
    Bootstrap(app)
    # EAM : Set limit on the number of items in cache (RAM)
    cache.init_app(app, config={'CACHE_TYPE': 'simple', 'CACHE_THRESHOLD': 1000})

    app.register_blueprint(frontend)

    app.json_encoder = MiniJSONEncoder

    nav.init_app(app)

    compress.init_app(app)
    htmlmin.init_app(app)

    return app


if __name__ == '__main__':
    create_app().start(debug=True)
