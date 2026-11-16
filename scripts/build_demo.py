#!/usr/bin/env python3
"""Inject results/*.json into demo/template.html -> demo/index.html (self-contained)."""
import json
import os

ROOT = os.path.join(os.path.dirname(__file__), "..")
load = lambda n: json.load(open(os.path.join(ROOT, "results", n + ".json")))
data = {"e1": load("e1_validation"), "e2": load("e2_error_dictionary"), "e3": load("e3_cost"), "e4": load("e4_noise")}
html = open(os.path.join(ROOT, "demo", "template.html"), encoding="utf-8").read()
assert "/*__DATA__*/null" in html
html = html.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))
open(os.path.join(ROOT, "demo", "index.html"), "w", encoding="utf-8").write(html)
print("  -> demo/index.html (%d kB)" % (len(html) // 1024))
