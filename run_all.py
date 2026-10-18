"""
run_all.py - run every step in order and stop at the first value that differs
from config.EXPECTED.

Outputs: numbers/key_numbers.json   snapshot of the main reported numbers
         numbers/env.txt            software-version sentence for the Methods

Usage:
    python run_all.py                 main analyses
    python run_all.py --biobank       also describe the serum-bank component (08)
    python run_all.py --symptoms      also run the respiratory-symptom description
"""
import sys
import json
import importlib.util
from datetime import datetime, timezone

import config as C

STEPS = [
    ("00_load_and_dedup", "raw records -> frozen participant files"),
    ("00b_extract_jihs", "JIHS 2023 survey, Table 3-1"),
    ("01_derive_and_validate", "flags from thresholds -> analysis file"),
    ("02_table1_characteristics", "Table 1"),
    ("03_analysis1_jihs_vs_pwoh", "Analysis 1: general population 2023 vs PWoH 2025"),
    ("04_analysis2_forest", "Analysis 2: multivariable logistic regression"),
    ("05_analysis3_mcnemar", "Analysis 3: PWH paired 2021 -> 2025"),
    ("06_cd4cd8_stratified", "PWH by baseline CD4/CD8 ratio"),
    ("07_table2_international", "Table 2: published comparators"),
    ("09_fig6_priming_generation", "PWoH by priming-vaccine generation"),
    ("10_hiv_history", "HIV disease history and 2025 seropositivity within PWH"),
]
OPTIONAL = {
    "--biobank": ("08_biobank_descriptive", "serum-bank component of PWoH set A"),
    "--symptoms": ("symptoms_descriptive", "respiratory symptoms at blood draw (descriptive)"),
}


def load(name):
    spec = importlib.util.spec_from_file_location(name, C.REPO / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def env_snapshot() -> dict:
    import platform
    import pandas as pd
    import numpy as np
    import scipy
    import matplotlib
    import statsmodels
    return {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "statsmodels": statsmodels.__version__,
    }


def main(argv):
    steps = list(STEPS) + [OPTIONAL[a] for a in OPTIONAL if a in argv]
    key = {"generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}

    for i, (name, desc) in enumerate(steps):
        print(f"\n\n{'#' * 72}\n# [{i + 1}/{len(steps)}] {name} - {desc}\n{'#' * 72}")
        try:
            mod = load(name)
            out = mod.main()
        except C.ExpectationError as e:
            print(f"\n\n{'!' * 72}")
            print(f"!! value differs from config.EXPECTED in {name}")
            print(f"!! {e}")
            print(f"{'!' * 72}")
            return 1
        except C.PrivateDataUnavailable as e:
            print(f"\n\n{'!' * 72}")
            print(f"!! {name} needs participant-level data that are not in this repository")
            print(f"!! {e}")
            print(f"{'!' * 72}")
            return 2
        if isinstance(out, dict):
            key.update(out)

    key["env"] = env_snapshot()
    key["frozen_n"] = {k: C.EXPECTED[k] for k in
                       ("n_pwh_2021", "n_pwh_2025", "n_pwoh_2025")}
    key["frozen_n"]["n_jihs_2023"] = C.EXPECTED["jihs_pdf_n"]
    key["figure_names"] = dict(C.FIGURES)
    key["match_algorithm"] = C.MATCH_ALGORITHM
    C.KEY_NUMBERS.write_text(json.dumps(key, ensure_ascii=False, indent=2))

    C.ENV_FILE.write_text(
        C.METHODS_STATS_SENTENCE.format(
            py=key["env"]["python"], pd=key["env"]["pandas"], np=key["env"]["numpy"],
            sp=key["env"]["scipy"], sm=key["env"]["statsmodels"],
            mpl=key["env"]["matplotlib"]) + "\n"
    )

    log = C.checks_log()
    print(f"\n\n{'=' * 72}")
    print(f"OK  {len(steps)} steps completed; {len(log)} checks against config.EXPECTED all passed")
    print(f"{'=' * 72}")
    print("  numbers/key_numbers.json   main reported numbers")
    print("  numbers/env.txt            software versions for the Methods")
    print("\n  next: Stata  do verify_in_stata.do")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
