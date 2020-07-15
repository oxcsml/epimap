"""Simple server so that we can visualise local files.

Run with `python3 run.py` and then browse to `localhost:8000`
"""

import http.server
import socketserver

PORT = 8000

Handler = http.server.SimpleHTTPRequestHandler

with socketserver.TCPServer(("", PORT), Handler) as httpd:
    print("serving on port", PORT)
    httpd.serve_forever()
