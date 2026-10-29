"""Sensitivity checks on frozen results; never overwrite the primary objects."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import h5py
from anndata.io import read_elem
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .publication_figure_state import project_paths, protein_paths, read_plot_object


def exact_difference(first, second):
    first, second = np.asarray(first, float), np.asarray(second, float)
    if min(len(first), len(second)) == 0 or not np.isfinite(np.r_[first, second]).all():
        raise ValueError("Both donor groups must be non-empty and finite.")
    pooled = np.r_[first, second]
    selected = np.asarray(list(itertools.combinations(range(len(pooled)), len(first))))
    sums = pooled[selected].sum(axis=1)
    null = sums / len(first) - (pooled.sum() - sums) / len(second)
    difference = float(first.mean() - second.mean())
    return difference, float((np.abs(null) >= abs(difference)-1e-12).mean()), len(null)


def profile_test(profiles, categories, seed=0, permutations=9999):
    """Permutation of independent donor labels on Hellinger profiles."""
    values = np.sqrt(profiles[categories].to_numpy(float))
    labels = profiles["condition"].astype(str).to_numpy()
    levels = np.unique(labels)
    if len(levels) < 2 or any((labels == level).sum() < 2 for level in levels):
        return {"n_units": len(labels), "r_squared": np.nan, "p_permutation": np.nan}
    total = float(np.sum((values-values.mean(axis=0))**2))
    def between(assignments):
        result = np.zeros(len(assignments))
        for level in levels:
            membership = assignments == level
            sums = membership.astype(float) @ values
            result += np.sum(sums*sums, axis=1) / membership.sum(axis=1)
        return result - len(values)*np.sum(values.mean(axis=0)**2)
    observed = between(labels[None, :])[0]
    rng = np.random.default_rng(seed)
    shuffled = np.stack([rng.permutation(labels) for _ in range(permutations)])
    null = between(shuffled)
    return {"n_units": len(labels), "r_squared": observed/total,
            "p_permutation": (1 + (null >= observed-1e-12).sum())/(permutations+1)}


def condition_tests(profiles, categories, seed=0, permutations=99999):
    labels = profiles.condition.astype(str).to_numpy()
    levels = np.unique(labels)
    rows = []
    for category in categories:
        # Reuse the primary notebook's deterministic seed and a longer prefix
        # of its permutation stream, avoiding arbitrary between-method seeds.
        digest = hashlib.sha256(f"{seed}|kruskal::primary::{category}".encode()).hexdigest()
        rng = np.random.default_rng(int(digest[:8], 16))
        shuffled = np.stack([rng.permutation(labels) for _ in range(permutations)])
        values = profiles[category].to_numpy(float)
        if np.allclose(values, values[0]):
            rows.append({"programme": category, "p_permutation": 1.0})
            continue
        ranks = stats.rankdata(values)
        _, counts = np.unique(values, return_counts=True)
        n = len(values)
        tie = 1 - np.sum(counts**3-counts)/(n**3-n)
        def statistic(assignments):
            term = sum(((assignments == level) @ ranks)**2 / (labels == level).sum() for level in levels)
            return (12*term/(n*(n+1))-3*(n+1))/tie
        observed = statistic(labels[None, :])[0]
        null = statistic(shuffled)
        rows.append({"programme": category,
                     "p_permutation": (1+(null >= observed-1e-12).sum())/(permutations+1)})
    result = pd.DataFrame(rows)
    result["p_holm"] = multipletests(result.p_permutation, method="holm")[1]
    result["n_permutations"] = permutations
    result["monte_carlo_se_unadjusted"] = np.sqrt(result.p_permutation*(1-result.p_permutation)/(permutations+1))
    return result


def donor_profiles(obs, calls, categories):
    frame = obs.loc[calls.index, ["sample", "patient", "condition", "analysis_unit"]].copy()
    frame["call"] = calls.astype(str)
    counts = frame.groupby("sample", observed=True).size()
    frame = frame[frame["sample"].isin(counts[counts >= 100].index)]
    sample = pd.crosstab([frame["sample"], frame.analysis_unit, frame.condition], frame.call, normalize="index")
    sample = sample.reindex(columns=categories, fill_value=0).reset_index()
    profiles = sample.groupby(["analysis_unit", "condition"], observed=True)[categories].mean().reset_index()
    if not np.allclose(profiles[categories].sum(axis=1), 1):
        raise ValueError("Incomplete programme composition.")
    return profiles


def doublet_sensitivity(epi, output):
    obs = epi.obs.copy()
    if obs.doublet_score.isna().any() or not obs.doublet_class.astype(str).eq("singlet").all():
        raise ValueError("This check expects retained singlets with finite doublet scores.")
    ror1 = epi[:, "ROR1"].layers["counts"].toarray().ravel()
    egfr = epi[:, "EGFR"].layers["counts"].toarray().ravel()
    summaries, donor_rows, cutoffs, programme_rows = [], [], [], []
    scored = ~obs.MP_programme_call.astype(str).isin(["not_scored", "normal_tissue"])
    categories = list(epi.uns["meta_programme_hallmark_order"]) + ["low_activity", "mixed_active_but_ambiguous"]
    for quantile in (1.0, .95, .90):
        limits = obs.groupby("sample", observed=True).doublet_score.quantile(quantile)
        keep = obs.doublet_score.le(obs["sample"].map(limits).astype(float)).to_numpy()
        for sample, limit in limits.items():
            mask = obs["sample"].eq(sample).to_numpy()
            cutoffs.append({"retained_score_quantile": quantile, "sample": sample, "score_cutoff": limit,
                            "n_before": int(mask.sum()), "n_after": int((mask & keep).sum())})
        for min_umi in (1, 2):
            frame = obs.loc[keep, ["condition", "patient"]].copy()
            frame["double_positive"] = ((ror1 >= min_umi) & (egfr >= min_umi))[keep]
            donors = frame.groupby(["condition", "patient"], observed=True).double_positive.agg(["size", "sum"]).reset_index()
            donors = donors[donors["size"] >= 100].copy()
            donors["prevalence"] = 100*donors["sum"]/donors["size"]
            donors["retained_score_quantile"], donors["min_umi"] = quantile, min_umi
            donor_rows.append(donors)
            tn = donors.loc[donors.condition.eq("TN"), "prevalence"]
            normal = donors.loc[donors.condition.eq("N"), "prevalence"]
            difference, p_value, n_perm = exact_difference(tn, normal)
            summaries.append({"retained_score_quantile": quantile, "min_umi": min_umi,
                              "tn_donors": len(tn), "normal_donors": len(normal),
                              "difference_pp": difference, "exact_permutation_p": p_value,
                              "n_permutations": n_perm})
        profiles = donor_profiles(obs, obs.loc[keep & scored, "MP_programme_call"], categories)
        programme_rows.append({"retained_score_quantile": quantile, **profile_test(profiles, categories)})
    pd.DataFrame(summaries).to_csv(output/"doublet_score_codetection_sensitivity.csv", index=False)
    pd.concat(donor_rows).to_csv(output/"doublet_score_donor_prevalence.csv", index=False)
    pd.DataFrame(cutoffs).to_csv(output/"doublet_score_cutoffs.csv", index=False)
    pd.DataFrame(programme_rows).to_csv(output/"doublet_score_programme_global.csv", index=False)
    return pd.DataFrame(summaries), pd.DataFrame(programme_rows)


def class_ranking_sensitivity(epi, tables, output):
    """Compare a coherent alternative aggregation; do not replace primary calls."""
    scores = pd.read_csv(tables/"meta_programme_class_scores.csv", index_col=0)
    activity = pd.read_csv(tables/"meta_programme_mp_activity_calls.csv", index_col=0)
    if not scores.index.equals(activity.index) or activity.isna().any().any():
        raise ValueError("Programme score/activity rows do not align.")
    groups = {"Cell Cycle": range(1, 5), "Stress or Hypoxia": range(5, 8), "Protein Regulation": range(8, 12),
              "Mesenchymal": range(12, 17), "Interferon and MHC II": range(17, 19), "Senescence": [19],
              "Oncogenic": [20], "Respiration": [21], "Secreted": [22, 23], "Cilia": [24],
              "Neural Lineage": range(25, 30), "Other Lineage": range(30, 33), "Haematological Lineage": range(33, 38)}
    ids = activity.columns.str.extract(r"^(MP\d+)", expand=False)
    mask = pd.DataFrame({category: activity.loc[:, ids.isin([f"MP{x}" for x in numbers])].any(axis=1)
                         for category, numbers in groups.items()}, index=scores.index)
    # Both the centring/scaling and ranking now use the same all-member statistic.
    ranked = scores.where(mask, -np.inf)
    n_active = mask.sum(axis=1)
    ordered = np.sort(ranked.to_numpy(float), axis=1)
    margin = np.full(len(scores), np.inf)
    multiple = n_active.ge(2).to_numpy()
    margin[multiple] = ordered[multiple, -1]-ordered[multiple, -2]
    alternative = ranked.idxmax(axis=1).where(n_active.gt(0), "low_activity")
    alternative = alternative.where(~(n_active.gt(0) & (margin < .25)), "mixed_active_but_ambiguous")
    primary = epi.obs.loc[scores.index, "MP_programme_call"].astype(str)
    categories = list(scores) + ["low_activity", "mixed_active_but_ambiguous"]
    pd.crosstab(primary.rename("primary"), alternative.rename("all_member_ranking")).to_csv(output/"class_ranking_confusion.csv")
    counts = pd.concat([primary.value_counts().rename("primary"), alternative.value_counts().rename("all_member_ranking")], axis=1).fillna(0)
    counts.to_csv(output/"class_ranking_assignment_counts.csv")
    pd.DataFrame({"primary": primary, "all_member_ranking": alternative}).to_csv(output/"class_ranking_cell_calls.csv")
    summaries, tests = [], []
    for name, calls in (("primary", primary), ("all_member_ranking", alternative)):
        profiles = donor_profiles(epi.obs, calls, categories)
        profiles.to_csv(output/f"class_ranking_donor_profiles_{name}.csv", index=False)
        summaries.append({"method": name, **profile_test(profiles, categories)})
        test = condition_tests(profiles, list(scores))
        test["method"] = name
        tests.append(test)
    summary = pd.DataFrame(summaries)
    summary["cell_assignment_agreement"] = float(primary.eq(alternative).mean())
    any_activity = primary.ne("low_activity")
    both_biological = primary.isin(scores.columns) & alternative.isin(scores.columns)
    summary["agreement_among_activity_supported"] = float(primary[any_activity].eq(alternative[any_activity]).mean())
    summary["agreement_among_both_biologically_assigned"] = float(primary[both_biological].eq(alternative[both_biological]).mean())
    summary.to_csv(output/"class_ranking_global_sensitivity.csv", index=False)
    pd.concat(tests).to_csv(output/"class_ranking_condition_tests.csv", index=False)
    return summary, counts


def correlation_sensitivity(frame, x_name, y_name, seed=0, permutations=99999, bootstrap=10000, subtype_adjust=False):
    complete = frame.dropna(subset=[x_name, y_name]).copy()
    x, y = (stats.rankdata(complete[name]) for name in (x_name, y_name))
    if subtype_adjust:
        design = np.c_[np.ones(len(complete)), pd.get_dummies(complete.subtype, drop_first=True).to_numpy(float)]
        x = x-design@np.linalg.lstsq(design, x, rcond=None)[0]
        y = y-design@np.linalg.lstsq(design, y, rcond=None)[0]
    x, y = x-x.mean(), y-y.mean()
    observed = float(np.dot(x, y)/(np.linalg.norm(x)*np.linalg.norm(y)))
    rng = np.random.default_rng(seed)
    groups = [np.flatnonzero(complete.subtype.to_numpy() == group) for group in complete.subtype.unique()] if subtype_adjust else [np.arange(len(y))]
    exceedances = 0
    for start in range(0, permutations, 1000):
        block = np.tile(y, (min(1000, permutations-start), 1))
        for indices in groups:
            for row in block:
                row[indices] = rng.permutation(y[indices])
        null = block@x/(np.linalg.norm(x)*np.linalg.norm(y))
        exceedances += int((np.abs(null) >= abs(observed)-1e-12).sum())
    result = {"n": len(complete), "rho": observed, "permutation_p": (exceedances+1)/(permutations+1)}
    if not subtype_adjust:
        raw_x, raw_y = complete[x_name].to_numpy(), complete[y_name].to_numpy()
        values = []
        for _ in range(bootstrap):
            indices = rng.integers(0, len(complete), len(complete))
            value = stats.spearmanr(raw_x[indices], raw_y[indices]).statistic
            if np.isfinite(value): values.append(value)
        result["ci_low"], result["ci_high"] = np.quantile(values, [.025, .975])
    return result


def protein_sensitivity(repository, output, project=None):
    protein_root = protein_paths(project)
    path = protein_root/"output/three_proteome_sample_scores.csv"
    scores = pd.read_csv(path)
    kosok = scores[scores.cohort.eq("Kosok")].copy()
    basal = kosok[kosok.subtype.isin(["Basal A", "Basal B"])].copy()
    if len(kosok) != 18 or len(basal) != 12:
        raise ValueError("The frozen cell-line panel changed.")
    pairs = [("ROR1–EGFR", "EGFR_z", "ROR1_z"), ("Pair–EMT", "pair_abundance_score", "EMT_score"),
             ("Pair–BCSC", "pair_abundance_score", "BCSC_score")]
    # Also rescore proteins within the basal subset: restricting observations
    # alone retains feature weights determined by the full 18-line panel.
    raw = pd.read_excel(protein_root/"data/mmc2.xlsx", sheet_name="S1A median log2 fold change", skiprows=4)
    prefix = "Normalized log2 fold-change: "
    columns = [column for column in raw if str(column).startswith(prefix)]
    matrix = raw[["Gene", *columns]].rename(columns={"Gene": "gene", **{c: str(c).removeprefix(prefix) for c in columns}})
    matrix = matrix.dropna(subset=["gene"])
    columns = list(matrix.columns[1:])
    matrix[columns] = matrix[columns].apply(pd.to_numeric, errors="coerce")
    matrix = matrix.groupby("gene")[columns].median().loc[:, basal.sample_id]
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    gmt = protein_root/"data/h.all.v2026.1.Hs.symbols.gmt"
    emt_genes = next(line.rstrip().split("\t")[2:] for line in gmt.read_text().splitlines()
                     if line.startswith("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION\t"))
    rescored = basal.set_index("sample_id").copy()
    for gene in ("ROR1", "EGFR"):
        rescored[f"{gene}_z"] = z.loc[gene]
    rescored["pair_abundance_score"] = rescored[["ROR1_z", "EGFR_z"]].mean(axis=1, skipna=False)
    for name, genes in (("EMT_score", emt_genes), ("BCSC_score", ["CD44", "ALDH1A1", "ALDH1A3", "ITGA6", "PROM1"])):
        values = z.loc[z.index.intersection(genes)]
        rescored[name] = values.mean().where(values.notna().sum() >= np.ceil(len(values)*.5))
    rescored.to_csv(output/"kosok_basal_rescored_samples.csv")
    rows = []
    for number, (association, x, y) in enumerate(pairs):
        rows.append({"analysis": "basal_subset_original_scores", "association": association,
                     **correlation_sensitivity(basal, x, y, seed=42+number)})
        rows.append({"analysis": "basal_subset_rescored", "association": association,
                     **correlation_sensitivity(rescored, x, y, seed=242+number)})
        rows.append({"analysis": "all_models_subtype_adjusted", "association": association,
                     **correlation_sensitivity(kosok, x, y, seed=142+number, subtype_adjust=True)})
    result = pd.DataFrame(rows)
    result["p_holm_sensitivity_family"] = multipletests(result.permutation_p, method="holm")[1]
    result.to_csv(output/"kosok_subtype_sensitivity.csv", index=False)
    missing = []
    for cohort, group in scores.groupby("cohort", sort=False):
        complete = group[["ROR1_z", "EGFR_z"]].notna().all(axis=1)
        for status, subset in (("paired_targets", group[complete]), ("missing_target", group[~complete])):
            content = pd.to_numeric(subset.get("tumour_content", pd.Series(dtype=float)), errors="coerce").dropna()
            missing.append({"cohort": cohort, "status": status, "n": len(subset),
                            "median_tumour_content": float(content.median()) if len(content) else np.nan})
    pd.DataFrame(missing).to_csv(output/"proteome_target_missingness.csv", index=False)
    return result


def barcode_audit(epi, project, repository, output):
    source = pd.read_csv(repository/"data/source_tumour_cluster_barcodes.csv.gz", dtype={"sample_code": str, "barcode": str})
    groups = {key: set(block.barcode) for key, block in source.groupby(["condition", "sample_code"])}
    rows = []
    for sample, obs in epi.obs.groupby("sample", observed=True):
        condition = str(obs.condition.iloc[0])
        match = __import__("re").search(r"(?:MH|AH|PM|SH|Tum)(\d+)", str(sample))
        code = match.group(1).zfill(4) if match else None
        key = (condition, code)
        barcodes = set(obs.index.str.split("-GSM").str[0])
        released = groups.get(key, set())
        head = project/"Adata_Saves/head"/f"{sample}.h5ad"
        retained = set()
        if head.exists() and released:
            with h5py.File(head) as handle:
                retained = set(read_elem(handle["obs"]).index.astype(str))
        overlaps = sorted(((len(barcodes & values), str(k)) for k, values in groups.items()), reverse=True)
        rows.append({"sample": sample, "condition": condition, "sample_code": code,
                     "source_group_available": bool(released), "n_source_barcodes": len(released),
                     "n_epithelial": len(barcodes), "source_in_epithelial": len(released & barcodes),
                     "source_in_retained_qc": len(released & retained),
                     "best_any_source_overlap": overlaps[0][0], "best_any_source_group": overlaps[0][1],
                     "excluded_from_transfer": str(sample)=="GSM4909296_ER-MH0001"})
    audit = pd.DataFrame(rows)
    audit["source_recovery_in_epithelium"] = audit.source_in_epithelial/audit.n_source_barcodes.replace(0, np.nan)
    audit.to_csv(output/"source_barcode_recovery_by_sample.csv", index=False)
    return audit


def run_checks(project=None, output_dir=None, sections=("doublet", "ranking", "protein", "barcodes")):
    project, repository = project_paths(project)
    output = Path(output_dir or repository/"results/publication_checks")
    output.mkdir(parents=True, exist_ok=True)
    object_path = project/"Adata_Saves/epithelial_antigen_coexpression.h5ad"
    epi = read_plot_object(object_path, genes=["ROR1", "EGFR"])
    results = {}
    if "doublet" in sections: results["doublet"] = doublet_sensitivity(epi, output)
    if "ranking" in sections: results["ranking"] = class_ranking_sensitivity(epi, project/"Dataframe_Saves", output)
    if "protein" in sections: results["protein"] = protein_sensitivity(repository, output, project=project)
    if "barcodes" in sections: results["barcodes"] = barcode_audit(epi, project, repository, output)
    provenance = {"object": object_path.name, "object_size_bytes": object_path.stat().st_size,
                  "object_mtime_ns": object_path.stat().st_mtime_ns, "sections": list(sections),
                  "primary_results_overwritten": False, "programme_calls_recomputed_after_filtering": False,
                  "notes": "Doublet filtering retains original programme calls; all added tests are sensitivities."}
    fingerprint = hashlib.sha256(pd.util.hash_pandas_object(epi.obs, index=True).to_numpy().tobytes())
    fingerprint.update(epi.layers["counts"].data.tobytes())
    fingerprint.update(epi.layers["counts"].indices.tobytes())
    fingerprint.update(epi.layers["counts"].indptr.tobytes())
    provenance["metadata_and_target_counts_sha256"] = fingerprint.hexdigest()
    provenance["check_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (output/"input_provenance.json").write_text(json.dumps(provenance, indent=2)+"\n")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--project")
    parser.add_argument("--output")
    parser.add_argument("--sections", nargs="+", default=["doublet", "ranking", "protein", "barcodes"])
    args = parser.parse_args()
    run_checks(args.project, args.output, tuple(args.sections))
