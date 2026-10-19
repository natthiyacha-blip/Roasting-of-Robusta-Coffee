#!/usr/bin/env python3
"""
HS-SPME-GC-MS curation pipeline (open-source, reproducible).
Inputs  (Agilent MassHunter Unknowns Analysis exports, one injection per condition):
  UnknownsReports/CEF/<COND>-AllBestHits.cef          -> RT, area, AMDIS spectrum, best NIST hit
  Mass spec/MFU-Sasi2025-BestHitDetail_*_<COND>.pdf    -> NIST match factor (MF) per component
  <COND>.D/tic_back.csv                                -> TIC trace (Fig. S1)
Steps
  1. pdftotext -layout on each Best-Hit PDF; parse ONLY the summary table (before the
     first 'Component RT:' detail page) -> RT, name, CAS, formula, MF, area.
  2. Parse CEF components (RT, area, deconvoluted spectrum) and attach MF by RT+name.
  3. Cross-sample alignment: greedy clustering of components by RT (|dRT| <= 0.06 min)
     AND spectral cosine (sqrt-intensity, m/z 29-330) >= 0.60, one component per sample
     per feature. Consensus annotation = highest-MF hit across samples.
  4. Curation: MF_max >= 80; remove chemically implausible hits (halogen/Si/B/P/D-bearing,
     silyl/borane/peroxy/nitro library artefacts, >C20); remove heavy hits (>=C13) that
     elute before 20 min on the polar column (elution-plausibility filter).
  5. Positional-isomer reassignment by elution order on WAX: 3-furanmethanol (RT 15.58)
     -> 2-furanmethanol; 3-furaldehyde (RT 11.30) -> furfural; C6H8N2 at RT 8.2-8.6 ->
     dimethylpyrazine isomer. Features with identical names within 0.06 min are merged.
  6. Chemical-class assignment by name/formula rules (see cls()).
  7. Retention-plausibility check: hits eluting before D-limonene with a nominal mass
     above 100 Da are relabelled 'Unidentified' (peak retained).
Usage: python gcms_curation.py <GCMS_root_folder>
Requires: poppler-utils (pdftotext), numpy, pandas, matplotlib
"""
import sys, os, re, glob, subprocess, xml.etree.ElementTree as ET
import numpy as np, pandas as pd
ROOT = sys.argv[1] if len(sys.argv) > 1 else "GCMS"
CONDS = ["HA-L", "HA-M", "HA-D", "MW-L", "MW-M", "MW-D"]
OUT = "gcms_outputs"; os.makedirs(OUT, exist_ok=True)

# ---- 1. PDF summary tables -> match factors
pat = re.compile(r'^\s*(\d+\.\d{4})\s+(.+?)\s+(\d{1,7}-\d{2}-\d)\s+(C[A-Za-z0-9]*)\s+(\d+\.\d)\s+([\d.]+)\s*$')
hits = []
for pdf in sorted(glob.glob(os.path.join(ROOT, "**", "*BestHitDetail*.pdf"), recursive=True)):
    cond = os.path.basename(pdf).split("_")[-1][:-4]
    if any(h["cond"] == cond for h in hits): continue
    txt = subprocess.run(["pdftotext", "-layout", pdf, "-"], capture_output=True, text=True).stdout.split("\n")
    end = next(i for i, l in enumerate(txt) if l.startswith("Component RT:")); last = None
    for l in txt[:end]:
        m = pat.match(l)
        if m:
            last = dict(cond=cond, rt=float(m[1]), name=m[2].strip(), formula=m[4], mf=float(m[5]), area=float(m[6])); hits.append(last); continue
        if not l.strip(): last = None; continue
        if last is not None and l.startswith(" " * 12) and not re.search(r"Page \d|Generated|Area|Component|Library Search", l):
            last["name"] = (last["name"] + ("" if last["name"].endswith("-") else " ") + l.strip()).strip()
hits = pd.DataFrame(hits); hits.to_csv(f"{OUT}/all_library_hits.csv", index=False)

# ---- 2. CEF components
comps = []
for f in sorted(glob.glob(os.path.join(ROOT, "**", "*-AllBestHits.cef"), recursive=True)):
    cond = os.path.basename(f)[:4]
    for c in ET.parse(f).getroot().iter("Compound"):
        loc, mol = c.find("Location"), c.find(".//Molecule")
        comps.append(dict(cond=cond, rt=float(loc.get("rt")), area=float(loc.get("a")),
                          name=mol.get("name") if mol is not None else None, formula=mol.get("formula") if mol is not None else None,
                          spec={round(float(p.get("x"))): float(p.get("y")) for p in c.iter("p")}))
cp = pd.DataFrame(comps)
def mf_of(r):
    s = hits[(hits.cond == r.cond) & ((hits.rt - r.rt).abs() < 0.003) & (hits.name == r["name"])]
    return s.mf.max() if len(s) else np.nan
cp["mf"] = cp.apply(mf_of, axis=1)
def vec(s):
    v = np.sqrt(np.array([s.get(m, 0.0) for m in range(29, 331)])); n = np.linalg.norm(v); return v / n if n else v
cp["v"] = cp.spec.map(vec); cp = cp.sort_values("rt").reset_index(drop=True)

# ---- 3. alignment
feats = []
for i, r in cp.iterrows():
    best, bs = None, 0
    for F in feats[-40:]:
        if abs(F["rt"] - r.rt) <= 0.06:
            s = float(F["v"] @ r.v)
            if s >= 0.60 and s > bs and r.cond not in F["m"]: best, bs = F, s
    if best is None: feats.append(dict(rt=r.rt, v=r.v.copy(), m={r.cond: i}))
    else: best["m"][r.cond] = i; k = len(best["m"]); best["rt"] = (best["rt"] * (k - 1) + r.rt) / k
rows = []
for F in feats:
    sub = cp.loc[list(F["m"].values())]; top = sub.sort_values("mf", ascending=False).iloc[0]
    row = dict(rt=F["rt"], name=top["name"], formula=top["formula"], mf_max=top.mf)
    for c in CONDS: row[c] = cp.loc[F["m"][c], "area"] if c in F["m"] else 0.0
    rows.append(row)
FT = pd.DataFrame(rows); FT.to_csv(f"{OUT}/feature_table_all.csv", index=False)

# ---- 4-5. curation
bad = re.compile(r'(F|Cl|Br|I|Si|B|Se|Sn|Ge|As|Hg|P|D)(?![a-z])')
def plaus(n, f):
    if not isinstance(f, str): return False
    if bad.search(re.sub(r'C|H|N|O|S', '', f)): return False
    if re.search(r'deuter|silyl|boran|fluor|chlor|brom|iod|anilino|trioxan|peroxy|nitro', n, re.I): return False
    return int(re.match(r'C(\d*)', f).group(1) or 1) <= 20
FT = FT[(FT.mf_max >= 80) & np.array([plaus(n, f) for n, f in zip(FT.name, FT.formula)])].copy()
FT["name"] = FT.name.str.replace(r"\s+", " ", regex=True)
nC = FT.formula.str.extract(r"^C(\d*)")[0].replace("", "1").astype(int)
FT = FT[~((nC >= 13) & (FT.rt < 20))].copy()
FT.loc[(FT.formula == "C6H8N2") & FT.rt.between(8.2, 8.6) & FT.name.str.contains("pyrimidine", case=False), "name"] = "Dimethylpyrazine isomer (C6H8N2)"
FT.loc[(FT.name == "3-Furanmethanol") & FT.rt.between(15.4, 15.8), "name"] = "2-Furanmethanol (furfuryl alcohol)"
FT.loc[(FT.name == "3-Furaldehyde") & FT.rt.between(11.1, 11.5), "name"] = "Furfural"
FT = FT.sort_values("rt"); g, grp, prev = 0, [], None
for _, r in FT.iterrows():
    if prev is not None and (r["name"] != prev["name"] or r.rt - prev.rt > 0.06): g += 1
    grp.append(g); prev = r
FT["g"] = grp
agg = {c: "sum" for c in CONDS}; agg.update(rt="mean", name="first", formula="first", mf_max="max")
FT = FT.groupby("g").agg(agg).reset_index(drop=True)

# ---- 6. classes
def cls(n, f):
    s = n.lower()
    if "caffeine" in s: return "Alkaloids"
    if "pyrazine" in s or "quinoxaline" in s: return "Pyrazines"
    if "pyridin" in s: return "Pyridines"
    if re.search(r"pyrrol|indole|indolizine", s): return "Pyrroles"
    if re.search(r"S\d*", re.sub(r"Si|Se|Sn", "", f)): return "Sulfur compounds"
    if re.search(r"phenol|guaiacol|catechol|cresol|methoxyphenyl|benzenediol", s): return "Phenols"
    if re.search(r"furan|furfur", s): return "Furans"
    if re.search(r"pyran|maltol|lactone|olide", s): return "Pyranones & lactones"
    if "acid" in s and "ester" not in s: return "Acids"
    if re.search(r"acetate|ester|oate\b", s): return "Esters"
    if re.search(r"al\b|aldehyde", s): return "Aldehydes"
    if re.search(r"one\b|dione|ketone", s): return "Ketones"
    if re.search(r"limonene|linalool|terpine|pinene|myrcene|ocimene", s): return "Terpenes"
    if re.search(r"ol\b|alcohol", s): return "Alcohols"
    if "N" in f: return "Other N-compounds"
    if re.fullmatch(r"C\d*H\d*", f): return "Hydrocarbons"
    return "Others"
FT["class"] = [cls(n, f) for n, f in zip(FT.name, FT.formula)]

# ---- 7. retention-plausibility relabelling
# On the polar WAX phase, compounds eluting before D-limonene (C10H16, 136 Da, RT 5.69 min)
# are restricted to small, highly volatile species. Library hits in this window with a
# nominal mass above 100 Da are inconsistent with their retention and are reported as
# "Unidentified" (peak kept for richness/diversity, excluded from named chemical classes).
AW = {"C": 12, "H": 1, "N": 14, "O": 16, "S": 32}
nom = lambda f: sum(AW.get(e, 0) * (int(k) if k else 1) for e, k in re.findall(r"([A-Z][a-z]?)(\d*)", f))
lim_rt = FT.loc[FT.name.str.contains("Limonene", case=False), "rt"].min()
unid = (FT.rt < lim_rt - 0.05) & (FT.formula.map(nom) > 100)
FT.loc[unid, "name"] = [f"Unidentified (RT {r:.2f} min)" for r in FT.loc[unid, "rt"]]
FT.loc[unid, "class"] = "Unidentified"; FT.loc[unid, "formula"] = "–"
print(f"Retention-plausibility check: {int(unid.sum())} features relabelled as Unidentified (RT < {lim_rt - 0.05:.2f} min, > 100 Da)")
FT[["rt", "name", "class", "formula", "mf_max"] + CONDS].to_csv(f"{OUT}/feature_table_curated.csv", index=False)
print(f"{len(hits)} library hits -> {len(FT)} curated aligned features; detected per sample:", (FT[CONDS] > 0).sum().to_dict())

# ---- Fig. S1  TIC overlay
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 8, "pdf.fonttype": 42})
# ---------------- Elsevier figure export ----------------
# Final widths: single column 90 mm, 1.5 column 140 mm, double column 190 mm.
# Each figure is resized so that its final (tight) width equals the target width,
# which keeps all text at its true point size (>= 7 pt), then exported as
# RGB TIFF (LZW, 1000 dpi), vector PDF and a 300 dpi PNG preview.
import io, logging
from PIL import Image
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
def save_elsevier(fig, name, width_mm, out=None, dpi=1000):
    out = out or OUT
    target = width_mm / 25.4
    for _ in range(12):
        fig.canvas.draw()
        w = fig.get_tightbbox(fig.canvas.get_renderer()).width + 0.04
        if 0 <= w - target < 0.002 * target: break
        target_adj = target * 1.0005
        s = target_adj / w
        fig.set_size_inches(fig.get_figwidth() * s, fig.get_figheight() * s)
        fig.tight_layout()
    kw = dict(bbox_inches="tight", pad_inches=0.02)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=dpi, **kw); buf.seek(0)
    im = Image.open(buf).convert("RGB")
    need = int(round(width_mm / 25.4 * dpi))                   # Elsevier minimum width in pixels
    if im.width < need:                                         # pad with white margin, no resampling
        canvas = Image.new("RGB", (need, im.height), "white"); canvas.paste(im, ((need - im.width) // 2, 0)); im = canvas
    im.save(f"{out}/{name}.tif", dpi=(dpi, dpi), compression="tiff_lzw")
    fig.savefig(f"{out}/{name}.pdf", **kw)
    fig.savefig(f"{out}/{name}.png", dpi=300, **kw)

COL = {"HA": "#eb6834", "MW": "#2a78d6"}
fig, axs = plt.subplots(3, 1, figsize=(7.2, 6), sharex=True)
for ax, d, lab in zip(axs, "LMD", ["Light", "Medium", "Dark"]):
    for r in ("HA", "MW"):
        f = glob.glob(os.path.join(ROOT, "**", f"{r}-{d}.D", "tic_back.csv"), recursive=True)
        if not f: continue
        t = pd.read_csv(f[0], skiprows=2, header=None, names=["t", "i"])
        ax.plot(t.t, t.i / 1e6 * (1 if r == "MW" else -1), color=COL[r], lw=0.7, label="Microwave (MW)" if r == "MW" else "Hot-air (HA, mirrored)")
    ax.axhline(0, color="#52514e", lw=0.5); ax.set_ylabel("TIC (×10$^6$)"); ax.set_title(lab, loc="left", fontsize=9, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{abs(v):g}"))
axs[0].legend(fontsize=7, frameon=False, loc="upper right"); axs[-1].set_xlabel("Retention time (min)"); axs[-1].set_xlim(2, 45)
fig.tight_layout(); save_elsevier(fig, "FigS1_TIC_mirror", 190); plt.close(fig)
