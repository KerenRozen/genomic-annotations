import json
import zlib
from pathlib import Path


def load_json(path: Path):
    return json.loads(
        zlib.decompress(
            path.read_bytes()
        ).decode('ascii'),
        object_hook=lambda d: {int(k): v for k, v in d.items()}
    )


def dump_json(path: Path, obj):
    path.write_bytes(
        zlib.compress(
            json.dumps(obj).encode('ascii')
        )
    )
