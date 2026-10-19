# Submission package — Innovative Food Science & Emerging Technologies (IFSET)

**Colour-Matched but Compositionally Different: An Iso-Chromatic Energy–Quality Comparison of a Rotary Microwave and a Hot-Air Roaster for Robusta Coffee**
Thanapon Saengsuwan, Chiangrai Rajabhat University

## Files and Editorial Manager item types

| Folder | File(s) | Upload as |
|---|---|---|
| 01_Manuscript | Manuscript_IFSET_MW_HA_Robusta.docx | Manuscript (title page, abstract, keywords, text, tables, figure captions, references; double-spaced, line numbers) |
| 02_Highlights | Highlights_IFSET.docx | Highlights (separate file; name contains "highlights") |
| 03_Figures | Figure_1 … Figure_6 (.tif, 1000 dpi; .pdf vector) | Figure (one file per figure) |
| 04_Tables | Tables_1-6.docx | Table (optional; tables are also in the manuscript) |
| 05_Supplementary_Material | Supplementary_Material.docx (Text S1, Fig. S1, Tables S1–S6), Figure_S1.tif, Table_S5_curated_volatile_features.xlsx | Supplementary material |
| 06_Code_and_Data | mw_roasting_analysis.py, gcms_curation.py, batch_level_data.csv, requirements.txt | Supplementary material (code and data) |

## Compliance with the IFSET Guide for Authors

| Requirement | Status |
|---|---|
| Research paper ≤ 8,000 words including allowance for figures and tables | Main text 4,960 words + abstract 234 words + 6 tables + 6 figures (≈ 7,500 words with 200 words per display item) |
| Abstract ≤ 250 words, abbreviations defined, no references | 234 words; MW, HA, TOST, TOPSIS defined |
| 1–7 keywords | 7 |
| Highlights: 3–5 bullets, ≤ 85 characters, separate file | 5 bullets, 76–80 characters |
| Title page: title, authors, affiliation with country, corresponding author e-mail | Included |
| Single-column .docx, double spacing, line numbers | Yes |
| SI units (equivalent SI given for other units) | kWh given with MJ equivalents (1 kWh = 3.6 MJ) |
| Equations editable, numbered consecutively, variables in italics | 10 Word (OMML) equations |
| Tables editable, cited, notes below, no vertical rules | Yes (three-line tables) |
| Figures as separate files, logical names, ≥ 1000 dpi for line art | Figure_1–6 TIFF at 1000 dpi (≥ 7,480 px full width; 5,512 px for 1.5 column) plus vector PDF with embedded fonts |
| Colour figures accessible to colour-vision deficiency | Orange/blue CVD-safe palette |
| Numbered sections (1, 1.1 …); abstract not numbered | Yes |
| CRediT statement, competing interests, funding, data statement | Included (ARDA funding statement) |
| GenAI declaration as a new section before the references | Included |
| Acknowledgements directly before the reference list | Yes |
| APA 7 references with DOIs, alphabetical | 74 references, all cited in text |
| Maximum 5–6 authors | 1 author listed |

### Figure files

| File | Pixels | Size (mm) |
|---|---|---|
| Figure_1.tif | 7488 × 5779 | 190 × 147 |
| Figure_2.tif | 7490 × 3224 | 190 × 82 |
| Figure_3.tif | 7480 × 2788 | 190 × 71 |
| Figure_4.tif | 7480 × 3537 | 190 × 90 |
| Figure_5.tif | 5512 × 6587 | 140 × 167 |
| Figure_6.tif | 7492 × 7409 | 190 × 188 |
| Figure_S1.tif | 7480 × 6203 | 190 × 158 |

## To complete before submission

1. **Funding:** ARDA statement included; add the grant/contract number if one exists.
2. **Research data (IFSET Option C):** deposit the data and code (folder 06 plus the raw GC-MS files) in a repository such as Mendeley Data, then add the DOI as a [dataset] reference and in the Data availability statement. Until then, the current statement explains how the data can be accessed.
3. **Declaration of interest** form generated with the Elsevier declarations tool (upload as .docx).
4. **Cover letter** (optional but recommended).

## Reproducing the analysis

```
pip install -r 06_Code_and_Data/requirements.txt
python 06_Code_and_Data/mw_roasting_analysis.py          # all statistics, Tables 1–6, Tables S1–S4, Figures 1–6
python 06_Code_and_Data/gcms_curation.py <GCMS_root>    # 259 curated features (Table S5), Figure S1
```
