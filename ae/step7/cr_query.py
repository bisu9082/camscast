"""Look up candidate references on Crossref by bibliographic query; print the best hits."""
import json, sys, time, urllib.parse, urllib.request
qs = [l.strip() for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
out = []
for q in qs:
    url = "https://api.crossref.org/works?rows=3&query.bibliographic=" + urllib.parse.quote(q)
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "refcheck/1.0"}), timeout=40) as r:
            items = json.load(r)["message"]["items"]
    except Exception as e:
        print("ERR", q, e); continue
    print("Q:", q)
    for it in items[:3]:
        t = (it.get("title") or [""])[0]
        a = (it.get("author") or [{}])[0].get("family", "")
        y = (it.get("issued", {}).get("date-parts") or [[None]])[0][0]
        c = (it.get("container-title") or [""])[0]
        print(f"   {it['DOI']} | {a} {y} | {t[:110]} | {c[:50]}")
    out.append(dict(q=q, hits=[dict(doi=i["DOI"], title=(i.get("title") or [""])[0]) for i in items[:3]]))
    time.sleep(0.6)
json.dump(out, open(sys.argv[1] + ".json", "w"), indent=1)
