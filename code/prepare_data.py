"""Export the three study workbooks as matrices without sample IDs or covariates."""

import argparse
import csv
from collections import Counter, defaultdict
from math import isfinite
from pathlib import Path

from openpyxl import load_workbook


DATA = Path(__file__).resolve().parents[1] / "data"
SAMPLES = {
    "CRC": {"C": 66, "H": 92, "P": 76},
    "METH": {"M": 170, "H": 50},
    "T2DM": {"D": 100, "H": 200},
}


def sheet_rows(workbook_path, sheet):
    workbook = load_workbook(workbook_path, read_only=True, data_only=True)
    try:
        return list(workbook[sheet].values)
    finally:
        workbook.close()


def meth_compounds(feature_names, feature_map, compound_names):
    if compound_names is None:
        with feature_map.open(encoding="utf-8-sig", newline="") as stream:
            return {row["Name"]: row["Compound"] for row in csv.DictReader(stream)
                    if row["Cohort"] == "METH"}

    synonyms = defaultdict(set)
    with compound_names.open(encoding="utf-8-sig") as stream:
        for line in stream:
            compound, aliases = line.rstrip("\r\n").split("\t", 1)
            for alias in aliases.split(";"):
                synonyms[alias.strip().casefold()].add(compound.removeprefix("cpd:"))
    mapping = {}
    for name in feature_names:
        matches = synonyms[name.casefold()]
        mapping[name] = next(iter(matches)) if len(matches) == 1 else ""
    return mapping


def make_crc(source):
    header, *records = sheet_rows(source / "CRC_data.xlsx", "Sheet1")
    samples = header[2:]
    features = [(str(compound).strip(), str(compound).strip(), name, abundances)
                for compound, name, *abundances in records]
    return [sample[0] for sample in samples], features


def make_meth(source, feature_map, compound_names):
    _, header, *records = sheet_rows(source / "MTBLS7807-Supplementary table.xlsx", "table2")
    samples = header[1:]
    name_to_compound = meth_compounds([r[0] for r in records], feature_map, compound_names)
    features = [(f"METH_{idx:03d}", name_to_compound[name], name, abundances)
                for idx, (name, *abundances) in enumerate(records, start=1)]
    return ["H" if sample.startswith("CTL") else "M" for sample in samples], features


def make_t2dm(source):
    header, *records = sheet_rows(
        source / "YLab_T2DM_FinalData_KEGG_pathways_metabolism_only.xlsx", "Data"
    )
    samples = header[3:]
    occurrences = Counter()
    features = []
    for name, compound, _, *abundances in records:
        occurrences[compound] += 1
        duplicate = occurrences[compound]
        feature = compound if duplicate == 1 else f"{compound}__{duplicate}"
        features.append((feature, compound, name, abundances))
    return [sample[0] for sample in samples], features


def validate(cohort, groups, features):
    feature_ids = {feature for feature, _, _, _ in features}
    expected_width = {"CRC": 113, "METH": 229, "T2DM": 146}
    expected_mapping = {"CRC": 113, "METH": 90, "T2DM": 144}
    compounds = {compound for _, compound, _, _ in features if compound}
    if (Counter(groups) != SAMPLES[cohort] or len(feature_ids) != len(features)
            or len(features) != expected_width[cohort] or len(compounds) != expected_mapping[cohort]):
        raise ValueError(f"{cohort} samples or feature mapping differ from the frozen inputs")
    for _, _, name, values in features:
        if len(values) != len(groups) or any(
            not isinstance(value, (int, float)) or not isfinite(value) or value <= 0
            for value in values
        ):
            raise ValueError(f"Invalid reported abundances for {cohort} feature {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Directory with the three study workbooks")
    parser.add_argument("--output", type=Path, default=DATA, help="Output directory")
    parser.add_argument("--feature-map", type=Path, default=DATA / "feature_map.csv",
                        help="Use the bundled frozen annotation map for METH")
    parser.add_argument("--kegg-compound-names", type=Path,
                        help="Optional frozen KEGG name/synonym TSV to rebuild the METH map")
    args = parser.parse_args()

    datasets = {
        "CRC": make_crc(args.source),
        "METH": make_meth(args.source, args.feature_map, args.kegg_compound_names),
        "T2DM": make_t2dm(args.source),
    }
    for cohort, (groups, features) in datasets.items():
        validate(cohort, groups, features)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "feature_map.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["Cohort", "Feature", "Compound", "Name"])
        for cohort, (groups, features) in datasets.items():
            ids, compounds, names, columns = zip(*features, strict=True)
            writer.writerows(zip([cohort] * len(ids), ids, compounds, names, strict=True))
            with (args.output / f"{cohort}.csv").open("w", encoding="utf-8", newline="") as matrix:
                rows = csv.writer(matrix)
                rows.writerow(["Group", *ids])
                rows.writerows(zip(groups, *columns, strict=True))


if __name__ == "__main__":
    main()
