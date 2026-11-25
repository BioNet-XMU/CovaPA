"""Run CovaPA, Global Test and Welch/Fisher ORA on the included cohort inputs."""

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))
from covapa_reference import covapa

CONTRASTS = {
    "CRC": (("C", "H"), ("C", "P"), ("P", "H")),
    "METH": (("M", "H"),),
    "T2DM": (("D", "H"),),
}


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, headings, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headings)
        writer.writerows(rows)


def cohort_inputs(name):
    data = ROOT / "data"
    paths = [data / name_file for name_file in (
        f"{name}.csv", "feature_map.csv", "pathway_memberships.csv"
    )]
    for path in paths:
        if not path.is_file():
            raise SystemExit(f"Missing data/{path.name}. See README.md for cohort inputs.")
    cohort = read_csv(paths[0])
    features = list(cohort[0])
    features.remove("Group")
    labels = np.array([row["Group"] for row in cohort])
    raw = np.array([[float(row[name]) for name in features] for row in cohort])
    if not np.isfinite(raw).all():
        raise ValueError(f"{name} contains nonfinite measurements.")
    mapping = {
        row["Feature"]: row["Compound"] for row in read_csv(paths[1])
        if row["Cohort"] == name
    }
    members = [row for row in read_csv(paths[2]) if row["Cohort"] == name]
    columns = defaultdict(list)
    for j, feature in enumerate(features):
        if compound := mapping[feature]:
            columns[compound].append(j)
    compounds = list(columns)
    lookup = {compound: j for j, compound in enumerate(compounds)}
    pathways = defaultdict(list)
    for row in members:
        pathways[row["KEGG_pathway"]].append(lookup[row["KEGG_compound"]])
    if name == "CRC":
        log_values = np.log(np.maximum(raw, 1e-12))
    else:
        # METH source Table 2 and the deposited T2DM matrix have positive
        # values for the mapped features used in CovaPA and Global Test.
        if any(np.any(raw[:, indices] <= 0) for indices in columns.values()):
            raise ValueError(f"{name} has nonpositive mapped abundances. Supply prepared inputs.")
        log_values = np.zeros_like(raw)
        np.log(raw, out=log_values, where=raw > 0)
    log_compounds = np.column_stack([log_values[:, idx].mean(axis=1) for idx in columns.values()])
    # Averaging duplicate records on the log scale becomes a geometric mean
    # on the reported-concentration scale used for the T2DM Global Test.
    gt = np.exp(log_compounds) if name == "T2DM" else log_compounds
    ora = raw if name == "T2DM" else log_values
    return labels, features, mapping, compounds, pathways, members, ora, log_compounds, gt


def compare(args, name, case, control, inputs):
    labels, features, mapping, compounds, pathways, members, ora, log_compounds, gt = inputs
    selected = np.isin(labels, [case, control])
    labels = labels[selected]
    if set(labels) != {case, control}:
        raise ValueError(f"Cannot find both groups for {name} {case}/{control}.")
    cov = covapa(log_compounds[selected], labels, pathways,
                 permutations=args.permutations, seed=args.seed)
    with tempfile.TemporaryDirectory(prefix=".run-", dir=ROOT) as folder:
        tmp = Path(folder)
        write_csv(tmp / "gt.csv", ["Group", *compounds], zip(labels, *gt[selected].T, strict=True))
        write_csv(tmp / "features.csv", ["Group", *features], zip(labels, *ora[selected].T, strict=True))
        write_csv(tmp / "map.csv", ["Feature", "Compound"], ((f, mapping[f]) for f in features))
        write_csv(tmp / "membership.csv", ["KEGG_pathway", "KEGG_compound"],
                  ((r["KEGG_pathway"], r["KEGG_compound"]) for r in members))
        # The METH archive reports a 91-compound background. The verified map
        # resolves 90 IDs, including all measured pathway members and observed hits.
        background = compounds if name == "METH" else sorted({r["KEGG_compound"] for r in members})
        universe_size = 91 if name == "METH" else len(background)
        (tmp / "background.txt").write_text("\n".join(background), encoding="ascii")
        inputs = [str((tmp / f).relative_to(ROOT)) for f in
                  ("gt.csv", "features.csv", "map.csv", "membership.csv", "background.txt")]
        command = [args.rscript, "code/run_comparators.R", *inputs,
                   str(args.permutations), str(args.seed), str(universe_size)]
        env = os.environ | ({"LC_ALL": "C"} if os.name == "nt" else {})
        try:
            result = subprocess.run(command, cwd=ROOT, env=env, check=True,
                                    stdout=subprocess.PIPE, text=True, encoding="utf-8")
        except FileNotFoundError as exc:
            raise SystemExit("Rscript is required. Use --rscript to give its full path.") from exc
        comparators = {r.pop("pathway"): r for r in csv.DictReader(result.stdout.splitlines())}
    dest = args.output / f"{name}_{case}_vs_{control}.csv"
    columns = ["Pathway", "Pathway_name", "Measured_members", "CovaPA_distance_squared",
               "CovaPA_p", "CovaPA_q", *next(iter(comparators.values()))]
    path_names = {r["KEGG_pathway"]: r["Pathway_name"] for r in members}
    write_csv(dest, columns, (
        [pid, path_names[pid], len(pathways[pid]), distance, p, q, *comparators[pid].values()]
        for pid, distance, p, q in cov
    ))
    print(f"{name} {case}/{control}: {len(labels)} samples, {len(pathways)} pathways -> {dest}")


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--cohort", choices=[*CONTRASTS, "all"], default="all")
    cli.add_argument("--permutations", type=int, default=999)
    cli.add_argument("--seed", type=int, default=20261119)
    cli.add_argument("--rscript", default="Rscript")
    cli.add_argument("--output", type=Path, default=ROOT / "results")
    args = cli.parse_args()
    if args.permutations < 1:
        cli.error("--permutations must be positive")
    args.output.mkdir(parents=True, exist_ok=True)
    for name in CONTRASTS if args.cohort == "all" else (args.cohort,):
        inputs = cohort_inputs(name)
        for case, control in CONTRASTS[name]:
            compare(args, name, case, control, inputs)


if __name__ == "__main__":
    main()
