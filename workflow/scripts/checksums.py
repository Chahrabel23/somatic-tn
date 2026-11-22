"""Write MD5 checksums of all inputs and a copy of the config used (traceability)."""
import hashlib

import yaml

sm = snakemake  # noqa: F821


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


with open(sm.output.md5, "w") as out:
    for f in dict.fromkeys(list(sm.input.fastq) + list(sm.input.ref)):
        out.write(f"{md5(f)}  {f}\n")
with open(sm.output.cfg, "w") as out:
    yaml.safe_dump(dict(sm.config), out, sort_keys=False)
