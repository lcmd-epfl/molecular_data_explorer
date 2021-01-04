# -*- coding: utf-8 -*-

# Run this app with `python mol_viewer.py` and
# visit http://127.0.0.1:8060/ in your web browser.

from app import create_app
import sys

data_file = sys.argv[1]
structures_dir = sys.argv[2]

app = create_app(data_file=data_file, structures_dir=structures_dir)

app.run_server(port=8060, debug=True)
