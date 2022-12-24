"""
1. Build index in memory and compare to the one on disk
2. Take a dictionary, dump and load then compare
"""
import hashlib
import json

from services.match_genes.build_index import build_index


def test__index_built_is_same():
    assert hashlib.sha1(
        json.dumps(
            build_index()
        ).encode("ascii")
    ).hexdigest() == "9946c88a59eadb3be45f88949ff639207cc14008"
