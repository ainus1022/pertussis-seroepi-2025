"""
00b_extract_jihs.py - general-population reference: JIHS fiscal-year 2023 survey,
Table 3-1 (age-specific distribution of anti-PT IgG).

The table is extracted from the published PDF report with `pdftotext -layout`
and written to data/jihs2023_pt_from_pdf.csv, which is shipped with this
repository. If the PDF is not available, the shipped CSV is validated instead.

Self-checks (either source):
  (a) all 71 rows (ages 0-69 and 70+): the titre bins sum to the row total
  (b) the column sums equal the Total row printed in the report
      [47, 329, 400, 624, 54, 18], n = 1472
"""
import re
import subprocess
import sys

import pandas as pd

import config as C

PDF_TOTAL_N = 1472
PDF_TOTAL_BINS = [47, 329, 400, 624, 54, 18]   # <1 / 1-4 / 5-9 / 10-49 / 50-99 / >=100
COLS = ["<1", "1-4", "5-9", "10-49", "50-99", ">=100"]

ROW_RE = re.compile(
    r"^\s*(70-|\d{1,2})\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"
    r"\s+([\d.]+)\s+([\d.]+)\s*$"
)


def validate(d: pd.DataFrame) -> None:
    if len(d) != 71:
        raise AssertionError(f"{len(d)} rows (expected 71: ages 0-69 and 70+)")
    bad = d[d[COLS].sum(axis=1) != d["Total"]]
    if len(bad):
        raise AssertionError(f"bins do not sum to the row total for ages {bad['age'].tolist()}")
    if int(d["Total"].sum()) != PDF_TOTAL_N:
        raise AssertionError(f"sum of totals {int(d['Total'].sum())} != {PDF_TOTAL_N}")
    if d[COLS].sum().tolist() != PDF_TOTAL_BINS:
        raise AssertionError(f"column sums {d[COLS].sum().tolist()} != {PDF_TOTAL_BINS}")


def extract() -> pd.DataFrame:
    txt = subprocess.run(["pdftotext", "-layout", str(C.SRC_JIHS_PDF), "-"],
                         capture_output=True, text=True, check=True).stdout

    # Tables 3-1 (PT) and 3-2 (FHA) have the same layout; the PT table is
    # identified by its Total row.
    anchor = re.compile(r"Total\s+1472\s+47\s+329\s+400\s+624\s+54\s+18\s")
    lines = txt.split("\n")
    starts = [i for i, ln in enumerate(lines) if anchor.search(ln)]
    if not starts:
        raise AssertionError("Total row of Table 3-1 (PT) not found; check the PDF layout")

    rows = {}
    for i in starts:
        for ln in lines[i:i + 90]:
            mm = ROW_RE.match(ln)
            if mm:
                g = mm.groups()
                rows.setdefault(g[0], dict(
                    age=g[0], Total=int(g[1]),
                    **{c: int(v) for c, v in zip(COLS, g[2:8])},
                    GMT=float(g[8]), GMT_log2=float(g[9])))
        if len(rows) >= 71:
            break

    d = pd.DataFrame(rows.values())
    d["age_num"] = pd.to_numeric(d["age"].str.replace("-", "", regex=False))
    d = d.sort_values("age_num").reset_index(drop=True)
    validate(d)
    d["ge10"] = d[["10-49", "50-99", ">=100"]].sum(axis=1)
    return d


def main():
    print("=" * 72)
    print("00b_extract_jihs - JIHS 2023 survey, Table 3-1 (anti-PT IgG)")
    print("=" * 72)

    if C.SRC_JIHS_PDF.exists():
        d = extract()
        d.to_csv(C.JIHS_FROM_PDF, index=False)
        print(f"  OK  extracted 71 rows from the PDF; row and column totals agree (n={PDF_TOTAL_N})")
    else:
        d = pd.read_csv(C.JIHS_FROM_PDF, dtype={"age": str})
        validate(d)
        print(f"  OK  PDF not available; shipped CSV validated (n={PDF_TOTAL_N})")

    adult = d[d["age_num"] >= 20]
    n, pos, high = int(adult["Total"].sum()), int(adult["ge10"].sum()), int(adult[">=100"].sum())
    print(f"\n  [age >= 20] n={n}  >=10: {pos} ({100 * pos / n:.1f}%)  "
          f">=100: {high} ({100 * high / n:.2f}%)")
    print(f"  -> {C.JIHS_FROM_PDF.relative_to(C.REPO)}")

    return {"jihs_pdf_n": n, "jihs_pdf_pos": pos, "jihs_pdf_high": high}


if __name__ == "__main__":
    main()
    sys.exit(0)
