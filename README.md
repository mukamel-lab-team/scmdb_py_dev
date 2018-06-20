# scmdb_py

## Development Setup
1. Clone the git repo.

  * `git clone https://github.com/mukamel-lab-team/scmdb_py_dev.git; cd scmdb_py_dev`

2. Set up virtual environment. 

  * `virtualenv venv -p python3`
  * `source venv/bin/activate`
  * `pip install --upgrade pip`
  * `pip install -r scmdb_py/requirements.txt`
  * `deactivate`

3. Set up default_config.py.
4. Run.
  * `bash run_dev.sh` *May need to change ports*
5. On your web browser, go to `localhost:<port>`
6. If running from a remote host (ie Banjo, Brainome), you must create an ssh tunnel  
from your local machine.
  * `ssh -NfL localhost:<port>:localhost:<port> <user_name>@<server_name>`


## Directory structure
```
scmdb_py_dev/
|-- scmdb_py.wsgi
|-- requirements.txt
|-- run_dev.sh
|-- scmdb_py/
|   |-- __init__.py
|   |-- frontend.py
|   |-- content.py
|   |-- assets.py
|   |-- default_config.py
|   |-- 
|   |-- 
```
