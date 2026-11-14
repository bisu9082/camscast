"""Write manuscript/references.bib holding only the entries cited in main_ae (from main_ae.aux)."""
import os
import re

MS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "manuscript")
aux = open(os.path.join(MS, "main_ae.aux")).read()
keys = []
for m in re.finditer(r"\\citation\{([^}]*)\}", aux):
    for k in m.group(1).split(","):
        if k.strip() and k.strip() not in keys:
            keys.append(k.strip())
bib = open(os.path.join(MS, "refs_ae.bib")).read()
entries = {}
for m in re.finditer(r"@(\w+)\{([^,]+),", bib):
    start = m.start()
    depth, i = 0, bib.index("{", start)
    for j in range(i, len(bib)):
        depth += {"{": 1, "}": -1}.get(bib[j], 0)
        if depth == 0:
            entries[m.group(2).strip()] = bib[start:j + 1]
            break
missing = [k for k in keys if k not in entries]
assert not missing, missing
open(os.path.join(MS, "references.bib"), "w").write("\n\n".join(entries[k] for k in keys) + "\n")
print(len(keys), "entries written")
