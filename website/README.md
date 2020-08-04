# Website

The website is static HTML, CSS and Javascript and technically can be run without a web-server. However, browsers block scripts from accessing local files, so we provide a very simple Python-based web-server to serve the local CSV files. 

- Run `python run.py`
- Navigate to localhost:8000 in your browser

If you want to view different Rt and CProj.csv files in the website directory, you can specify these as query parameters:
`http://localhost:8000/?rt=myrt.csv&cproj=mycproj.csv`

If you need to change the port, update the constant in run.py.
