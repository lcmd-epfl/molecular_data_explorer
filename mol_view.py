# -*- coding: utf-8 -*-

# Run this app with `python mol_viewer.py` and
# visit http://127.0.0.1:8060/ in your web browser.

from app import create_app
import sys

app = create_app()

app.run_server(port=8060, debug=True)
