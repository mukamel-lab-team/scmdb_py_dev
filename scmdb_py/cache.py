"""Flask-Cache instance."""
from flask_caching import Cache

cache = Cache(config={'CACHE_TYPE': 'simple', 'CACHE_THRESHOLD': 1000})
