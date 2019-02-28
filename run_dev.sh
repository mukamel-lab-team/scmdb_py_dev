source venv/bin/activate

# Compile the dissections html
php -f scmdb_py/templates/components/CEMBA_dissections.php > scmdb_py/templates/components/CEMBA_dissections.html

# Serve web app using flask on port 50xx
flask --app=scmdb_py dev -p 5033
