source venv/bin/activate

# Compile the dissections html
php -f scmdb_py/templates/components/CEMBA_dissections.php > scmdb_py/templates/components/CEMBA_dissections.html

# Serve web app using flask on port 50xx
#export FLASK_APP=scmdb_py/frontend.py
#export FLASK_ENV=development
#export FLASK_RUN_PORT=5033
#flask run

flask --app=scmdb_py dev -p 5033
