# scmdb_py

## Development Setup
1. Clone the git repo.

   * `https://github.com/JoeyHou/scmdb_py_dev_team_copy.git; cd scmdb_py_dev`

2. Set up virtual environment.  Note that scmdb_py currently works with python3.5

   * `virtualenv venv -p python3.5`
   * `source venv/bin/activate`
   * `pip install --upgrade pip`
   * `pip install -r requirements.txt`
   * `deactivate`

3. Create directory with proper permissions for user-login.sqlite.
   * `mkdir scmdb_py/tmp`
   * `chmod 777 scmdb_py/tmp`
   * Find and move user-login.sqlite file to scmdb_py/tmp and change permissions to read+write+execute on all.
4. Set up default_config.py.
   * Leave `APPLICATION_ROOT` blank
   * *Make sure to not push default_config.py to github unless you have removed  
all sensitive information!*
5. If the machine you are running this on does NOT have mysql installed, uncomment the lines `import pymsysql` and `pymysql.install_as_MySQLdb()` in \__init__.py.
6. Run.
   * `bash run_dev.sh` *May need to change ports*
7. On your web browser, go to `localhost:<port>`
8. If running from a remote host (ie Banjo, Brainome), you must create an ssh tunnel  
from your local machine.
   * `ssh -NfL localhost:<port>:localhost:<port> <user_name>@<server_name>`

## Deployment (requires sudo access)
1. Follow steps 1-5 of development setup.
2. Move whole directory to `/var/www/`
   * `cd ..`
   * `sudo mv scmdb_py_dev /var/www/scmdb_py`
3. Set `APPLICATION_ROOT` in default_config.py.
   * This sets the url root for the portal.
   * ex. To host the portal on website.com/portal, set `APPLICATION_ROOT='/portal'`.
4. Edit .wsgi file to point to proper file paths.
5. Edit .conf file in `/etc/apache2/sites-enabled` to point to proper file paths
   * 000-default.conf (http configuration) and default-ssl.conf (https configuration)
6. Make sure you have no syntax errors in the apache2 conf.
   * `apachectl configtest`
7. Make sure the proper version of mod_wsgi is installed (must be python3)
   * `ldd /usr/lib/apache2/modules/mod_wsgi.so`
   * If not installed or compiled with the wrong version of python, first remove the old package using apt,   
then install `sudo apt install libapache2-mod-wsgi-py3`
8. Restart apache2
   * `sudo service apache2 restart`


## Troubleshooting deployment setup
1. Read the error log
   * `sudo less /var/log/apache2/brainome-error_log`
2. Double check permissions for the /static and /tmp directory.
3. Make sure you did steps 3-6 of deployment setup.
   * mod_wsgi only works if it has been compiled with the same version of python as the virtualenv.
4. Check if default_config.py is properly setup.



## Directory structure
```
scmdb_py_dev/
|-- scmdb_py.wsgi                           *WSGI script file for hosting via apache *Don't Touch*
|-- requirements.txt                        *list of required python packages
|-- run_dev.sh
|-- scmdb_py/
|   |-- __init__.py                         *Application factory (setup)
|   |-- frontend.py                         *responsible for all views (handles URL requests)
|   |-- content.py                          *all server side data querying and plot generation
|   |-- assets.py                           *gathers all javascript files in assets directory
|   |-- default_config.py                   *Configuration file for Flask. (info for MySQL, email, etc.)
|   |-- assets/                             *All your .js and .css files go here
|   |   |-- scripts/                        *Individual javascript files
|   |   |   |-- app.js                      
|   |   |   |-- customview.js               *Frontend functions for data browser (plotting, sending queries)
|   |   |   |-- request_new_ensemble.js     *Frontend functions for request new ensemble page
|   |   |   |-- tabular_dataset_rs1.js      *Frontend functions for summary (metadata) pages
|   |   |   |-- tabular_dataset_rs2.js
|   |   |   |-- tabular_ensemble.js
|   |   |   |-- vendor/                     *All downloaded javascript libraries
|   |   |   |   |-- ...
|   |   |-- styles/                         *Individual css files
|   |   |   |-- app_base.css                *Style for navigation bar.
|   |   |   |-- browser.css                 *Style for data browser
|   |   |   |-- vendor/                     *All downloaded css files
|   |   |   |   |-- ...
|   |-- templates/                          *HTML files 
|   |   |-- base.html                       *Base template file for ensembleview.html
|   |   |-- ensembleview.html               *Main HTML file for data browser page of portal
|   |   |-- request_new_ensemble.html
|   |   |-- tabular_dataset_rs1.html
|   |   |-- tabular_dataset_rs2.html
|   |   |-- tabular_ensemble.html
|   |   |-- accounts/                       *HTML files for account/login pages
|   |   |-- components/                     *Individual componenents of ensembleview.html
|   |   |   |-- search_options.html         *Gene methylation search box
|   |   |   |-- mch_box.html                *Methylation box plot / heatmap
|   |   |   |-- mch_correlated_genes.html   *Methylation box plot / heatmap
|   |   |   |-- mch_scatter.html            *Methylation cluster plot
|   |   |-- partials/                       *HTML code to be included in all pages
|   |   |   |-- _flashes.html               *flash message settings
|   |   |   |-- _head.html                  *included in <head> of all pages
|   |   |-- macros/                         *Macros (functions) used by jinja to generate HTML
|   |   |   |-- nav_macros.html             *Edit this to change navigation bar items
|   |-- static/                             *All static files (.js files generated by Flask, images) *Don't Touch*
|   |-- tmp/
|   |   |-- user-login.sqlite               *User login information. DON'T PUSH TO GITHUB
```

## How requests are handled
1. User clicks a button on page. (ie. Search for the gene Sox6)
2. Javascript event handler detects an event and runs the appropriate function that  
sends an AJAX request via HTML GET or POST method.
3. frontend.py handles the GET/POST method and runs a function from content.py.
4. content.py queries MySQL, shapes the data using pandas and numpy, generates  
HTML for plots, and returns data as a list or dict to frontend.py.
5. frontend.py sends the data back to the client in JSON format.
6. Javascript updates what the user sees on screen. 

