source /srv/venv/bin/activate
cd /srv/
echo `date` Entered virtual environment.
echo Flask at: `which flask`

flask --app=scmdb_py serve --port 5000

echo `date` Killed / quit.
deactivate
