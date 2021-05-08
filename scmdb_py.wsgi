activate_this = '/var/www/scmdb_py/venv/bin/activate_this.py'
with open(activate_this) as file_:
    exec(file_.read(), dict(__file__=activate_this))

import os
import sys
import logging

logging.basicConfig(stream=sys.stderr)

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path)

#from scmdb_py import app as application
from scmdb_py import create_app

config_path = os.path.join(dir_path, 'scmdb_py/default_config.py')
application = create_app(config_path)
application.secret_key = 'cndd_ddnc'
