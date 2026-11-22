"""Step 7, layers 1-3 of the citation check for manuscript/references.bib.

Layer 1 arXiv (entries whose DOI is an arXiv DOI), layer 2 Crossref (DOI resolves; title, first
author and year agree with the bibliography entry), layer 3 Semantic Scholar (the work is
indexed under the same DOI and title). Writes step7/citation_layers.json.
Public endpoints only; no personal data is sent. A Semantic Scholar key is read from S2_KEY if set.
"""
import difflib
import json
import os
import re
import time
import unicodedata
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BIB = os.path.join(HERE, "..", "manuscript", "references.bib")
OUT = os.path.join(HERE, "..", "step7", "citation_layers.json")
UA = {"User-Agent": "citation-check/1.0"}


def get(url, tries=4):
    err = ""
    h = dict(UA)
    if "semanticscholar" in url and os.environ.get("S2_KEY"):
        h["x-api-key"] = os.environ["S2_KEY"]          # key from the environment only, never stored
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=40) as r:
                return r.read().decode("utf-8")
        except Exception as e:  # noqa: BLE001
            err = str(e)
            time.sleep(3 * (k + 1))
    return "ERROR: " + err


def entries():
    """Split the file at each '@' entry so that one malformed entry cannot swallow the next."""
    bib = open(BIB).read()
    out = []
    for chunk in re.split(r"\n(?=@\w+\{)", "\n" + bib):
        m = re.match(r"\s*@(\w+)\{([^,]+),(.*)", chunk, re.S)
        if not m:
            continue
        f = {}
        for fm in re.finditer(r"(\w+)\s*=\s*\{((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)\}", m.group(3)):
            f[fm.group(1).lower()] = " ".join(fm.group(2).split())
        f["key"], f["type"] = m.group(2).strip(), m.group(1).lower()
        out.append(f)
    return out


def norm(s):
    s = re.sub(r"\\[a-zA-Z]+|[{}$_^\\]", "", s or "")
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def sim(a, b):
    return round(difflib.SequenceMatcher(None, norm(a), norm(b)).ratio(), 3)


def first_family(author_field):
    a = (author_field or "").split(" and ")[0]
    fam = a.split(",")[0] if "," in a else (a.split()[-1] if a.split() else "")
    return norm(fam)


def main():
    only = [k for k in os.environ.get("ONLY", "").split(",") if k]
    prev = json.load(open(OUT)) if only and os.path.exists(OUT) else []
    res = [r for r in prev if r["key"] not in only]
    for e in entries():
        if only and e["key"] not in only:
            continue
        doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", (e.get("doi") or "").strip())
        r = dict(key=e["key"], title=e.get("title"), year=e.get("year"), doi=doi,
                 first_author=first_family(e.get("author")))
        m = re.search(r"arxiv\.(\d{4}\.\d{4,5})", doi, re.I)
        if m:
            x = get("https://export.arxiv.org/api/query?id_list=" + m.group(1))
            t = re.search(r"<entry>.*?<title>(.*?)</title>", x, re.S)
            r["arxiv"] = dict(id=m.group(1), found=bool(t),
                              title_sim=sim(t.group(1), e.get("title", "")) if t else None)
        if not doi:
            r["crossref"] = dict(found=False, error="no DOI")
            res.append(r)
            print(f"{r['key']:22s} NO DOI", flush=True)
            continue
        x = get("https://api.crossref.org/works/" + urllib.parse.quote(doi))
        if x.startswith("ERROR"):
            # DataCite DOIs (e.g. arXiv) are not in Crossref: check that the DOI resolves
            y = get("https://api.datacite.org/dois/" + urllib.parse.quote(doi))
            if not y.startswith("ERROR"):
                a = json.loads(y)["data"]["attributes"]
                t = (a.get("titles") or [{}])[0].get("title", "")
                r["crossref"] = dict(found=True, registry="DataCite", title=t,
                                     title_sim=sim(t, e.get("title", "")),
                                     years=[a.get("publicationYear")],
                                     year_match=str(a.get("publicationYear")) == str(e.get("year")),
                                     first_author_match=r["first_author"] in norm(
                                         json.dumps(a.get("creators", [])[:1])))
            else:
                r["crossref"] = dict(found=False, error=x[:120])
        else:
            msg = json.loads(x)["message"]
            ct = (msg.get("title") or [""])[0]
            fa = norm((msg.get("author") or [{}])[0].get("family", ""))
            yrs = sorted({msg[k]["date-parts"][0][0] for k in
                          ("published-print", "published-online", "issued", "published")
                          if msg.get(k, {}).get("date-parts") and msg[k]["date-parts"][0][0]})
            r["crossref"] = dict(found=True, registry="Crossref", title=ct,
                                 title_sim=sim(ct, e.get("title", "")),
                                 first_author=fa,
                                 first_author_match=bool(fa) and (fa == r["first_author"]
                                                                  or fa in r["first_author"]
                                                                  or r["first_author"] in fa),
                                 years=yrs, year_match=str(e.get("year")) in {str(y) for y in yrs},
                                 container=(msg.get("container-title") or [""])[0])
        time.sleep(0.8)
        y = get("https://api.semanticscholar.org/graph/v1/paper/DOI:" + urllib.parse.quote(doi)
                + "?fields=title,year,abstract")
        if y.startswith("ERROR"):
            r["s2"] = dict(found=False, error=y[:80])
        else:
            j = json.loads(y)
            r["s2"] = dict(found=True, title_sim=sim(j.get("title", ""), e.get("title", "")),
                           year=j.get("year"), has_abstract=bool(j.get("abstract")))
        time.sleep(1.6)
        res.append(r)
        cr, s2 = r["crossref"], r["s2"]
        print(f"{r['key']:22s} {cr.get('registry','-'):8s} tsim={cr.get('title_sim')} "
              f"fa={cr.get('first_author_match')} yr={cr.get('year_match')} "
              f"s2={s2.get('found')} s2sim={s2.get('title_sim')}", flush=True)
    order = [e["key"] for e in entries()]
    res.sort(key=lambda r: order.index(r["key"]) if r["key"] in order else len(order))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
