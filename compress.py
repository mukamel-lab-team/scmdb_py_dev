"""Instances of modules used to compress transfers."""
from flask_compress import Compress
from flask.ext.htmlmin import HTMLMIN

compress = Compress()
htmlmin = HTMLMIN()
