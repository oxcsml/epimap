# Website

The website is static HTML, CSS and Javascript and technically can be run without a web-server. However, browsers block scripts from accessing local files, so we provide a very simple Python-based web-server to serve the local CSV files. 

- Run `python run.py`
- Navigate to localhost:8000 in your browser

If you need to change the port, update the constant in run.py.
