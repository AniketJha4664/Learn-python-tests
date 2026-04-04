import csv
import math
import random
from dataclasses import dataclass

RANDOM_SEED = 42
random.seed(RANDOM_SEED)


@dataclass
class GeneDiscoveryResult:
    accuracy: float
    auc_approx: float
    top_driver_genes: list
    mutation_clusters: list


def sigmoid(x):
    return 1.0 / (1.0 + math.exp(-x))


def mean(values):
    return sum(values) / len(values) if values else 0.0


def stdev(values):
    m = mean(values)
    var = sum((x - m) ** 2 for x in values) / max(1, len(values) - 1)
    return math.sqrt(var)


def point_biserial(feature, labels):
    pos = [f for f, y in zip(feature, labels) if y == 1]
    neg = [f for f, y in zip(feature, labels) if y == 0]
    n1, n0 = len(pos), len(neg)
    if n1 == 0 or n0 == 0:
        return 0.0
    s = stdev(feature)
    if s == 0:
        return 0.0
    return ((mean(pos) - mean(neg)) / s) * math.sqrt((n1 * n0) / ((n1 + n0) ** 2))


def roc_auc_approx(scores, labels):
    paired = sorted(zip(scores, labels), key=lambda x: x[0])
    pos_ranks = 0
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    for i, (_, y) in enumerate(paired, start=1):
        if y == 1:
            pos_ranks += i
    if n_pos == 0 or n_neg == 0:
        return 0.5
    return (pos_ranks - (n_pos * (n_pos + 1) / 2)) / (n_pos * n_neg)


def generate_patient_genomic_data(n_patients=400, n_genes=80):
    genes = [f"GENE_{i:03d}" for i in range(n_genes)]
    driver_idx = [3, 7, 16, 25, 41, 58, 73]

    expressions = []
    labels = []
    mutations = []

    for _ in range(n_patients):
        cancer_signal = random.gauss(0, 1)
        y = 1 if sigmoid(1.3 * cancer_signal) > 0.5 else 0
        sample = [random.gauss(0, 1) for _ in range(n_genes)]

        for i in driver_idx[:4]:
            sample[i] += 1.4 * cancer_signal
        for i in driver_idx[4:]:
            sample[i] -= 1.1 * cancer_signal

        mut = [1 if random.random() < 0.05 else 0 for _ in range(n_genes)]
        for i in driver_idx:
            mut[i] = 1 if random.random() < 0.35 else 0

        expressions.append(sample)
        labels.append(y)
        mutations.append(mut)

    return genes, expressions, labels, mutations, [genes[i] for i in driver_idx]


def stage1_gene_detection(genes, expressions, labels, mutations):
    # Rank genes by association to tumor/normal labels
    importances = []
    for j, gene in enumerate(genes):
        feature = [row[j] for row in expressions]
        score = abs(point_biserial(feature, labels))
        importances.append((gene, score))
    importances.sort(key=lambda x: x[1], reverse=True)
    top_genes = importances[:12]

    # Simple weighted classifier using top genes
    selected = [genes.index(g) for g, _ in top_genes[:8]]
    weights = [point_biserial([row[i] for row in expressions], labels) for i in selected]

    probs = []
    preds = []
    for row in expressions:
        z = sum(w * row[i] for w, i in zip(weights, selected))
        p = sigmoid(z)
        probs.append(p)
        preds.append(1 if p >= 0.5 else 0)

    accuracy = sum(1 for p, y in zip(preds, labels) if p == y) / len(labels)
    auc = roc_auc_approx(probs, labels)

    # K-means on mutation matrix (binary features)
    k = 3
    centroids = [mutations[i][:] for i in random.sample(range(len(mutations)), k)]
    assignments = [0] * len(mutations)

    for _ in range(8):
        for idx, row in enumerate(mutations):
            dists = [sum((a - b) ** 2 for a, b in zip(row, c)) for c in centroids]
            assignments[idx] = dists.index(min(dists))
        for c in range(k):
            members = [mutations[i] for i, a in enumerate(assignments) if a == c]
            if members:
                centroids[c] = [1 if mean([m[j] for m in members]) >= 0.5 else 0 for j in range(len(genes))]

    cluster_counts = []
    for c in range(k):
        count = sum(1 for a in assignments if a == c)
        cluster_counts.append((c, count))

    return GeneDiscoveryResult(
        accuracy=accuracy,
        auc_approx=auc,
        top_driver_genes=top_genes,
        mutation_clusters=cluster_counts,
    )


def stage2_3_drug_design():
    # Candidate library with descriptor-like features
    candidates = []
    for i in range(180):
        mw = random.gauss(415, 65)
        lip = random.gauss(3.5, 0.9)
        hbond = random.randint(1, 10)
        docking = random.gauss(-8.2, 1.3)
        pathway = random.uniform(0.25, 1.0)

        efficacy = (
            72
            + 7 * pathway
            - 2.3 * abs((mw - 420) / 70)
            - 5 * (lip - 3.2) ** 2
            - 0.8 * abs(hbond - 5)
            - 2.5 * max(docking + 8, 0)
            + random.gauss(0, 1.6)
        )
        tox_prob = sigmoid(-2.6 + 0.014 * (mw - 450) + 0.58 * (lip - 4) + 0.2 * (hbond - 6))
        score = efficacy - 28 * tox_prob

        candidates.append({
            "id": f"CAND_{i+1:03d}",
            "mw": mw,
            "lip": lip,
            "hbond": hbond,
            "docking": docking,
            "pathway": pathway,
            "efficacy": efficacy,
            "tox_prob": tox_prob,
            "score": score,
        })

    candidates.sort(key=lambda x: x["score"], reverse=True)
    top = candidates[:8]

    # Generative optimization: perturb top candidates and keep good designs
    generated = []
    for parent in top:
        for _ in range(35):
            mw = min(700, max(180, parent["mw"] + random.gauss(0, 18)))
            lip = min(7, max(0.3, parent["lip"] + random.gauss(0, 0.25)))
            hbond = int(min(12, max(1, round(parent["hbond"] + random.gauss(0, 1)))))
            docking = parent["docking"] + random.gauss(0, 0.35)
            pathway = min(1, max(0, parent["pathway"] + random.gauss(0, 0.05)))

            efficacy = 72 + 7 * pathway - 2.3 * abs((mw - 420) / 70) - 5 * (lip - 3.2) ** 2 - 0.8 * abs(hbond - 5) - 2.5 * max(docking + 8, 0)
            tox_prob = sigmoid(-2.6 + 0.014 * (mw - 450) + 0.58 * (lip - 4) + 0.2 * (hbond - 6))

            if efficacy >= 74 and tox_prob <= 0.35:
                generated.append({
                    "id": f"GEN_{len(generated)+1:03d}",
                    "efficacy": efficacy,
                    "tox_prob": tox_prob,
                    "docking": docking,
                    "pathway": pathway,
                })

    generated.sort(key=lambda x: (x["efficacy"], -x["tox_prob"]), reverse=True)
    return top, generated[:10]


def stage4_delivery(n_patients=12):
    plans = []
    for i in range(n_patients):
        age = random.randint(35, 82)
        tumor = random.uniform(0.2, 0.95)
        metastatic = 1 if random.random() < 0.45 else 0
        kidney = random.uniform(45, 115)
        immune = random.uniform(0.1, 0.95)

        dose = 65 + 33 * tumor + 12 * metastatic - 0.2 * age - 0.17 * (100 - kidney) + 8 * immune
        dose = max(35, min(130, dose))

        if metastatic == 1 and tumor > 0.65:
            delivery = "Nanoparticle-targeted IV"
        elif immune > 0.65:
            delivery = "Antibody-drug conjugate"
        else:
            delivery = "Localized targeted therapy"

        plans.append({
            "patient_id": f"PT_{i+1:03d}",
            "tumor_burden": tumor,
            "metastatic": metastatic,
            "immune_score": immune,
            "recommended_dose_mg": dose,
            "delivery_mechanism": delivery,
        })
    return plans


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    print("=== Oncology AI Pipeline Demo (No external dependencies) ===")

    genes, expressions, labels, mutations, truth = generate_patient_genomic_data()
    stage1 = stage1_gene_detection(genes, expressions, labels, mutations)

    print("\n[1] Gene Detection & Target Identification")
    print(f"Classification accuracy: {stage1.accuracy:.3f}")
    print(f"ROC-AUC (approx): {stage1.auc_approx:.3f}")
    print("Top inferred genes:")
    for gene, imp in stage1.top_driver_genes:
        print(f"  {gene:>10} | importance={imp:.4f}")
    print("Mutation clusters:")
    for c, n in stage1.mutation_clusters:
        print(f"  Cluster {c}: {n} patients")
    print("Synthetic ground-truth drivers:", ", ".join(truth))

    top_candidates, generated = stage2_3_drug_design()
    print("\n[2] Oncogene Targeting Strategy")
    print("Top optimized candidates:")
    for c in top_candidates:
        print(
            f"  {c['id']} | efficacy={c['efficacy']:.2f} | tox={c['tox_prob']:.3f} "
            f"| docking={c['docking']:.2f} | pathway={c['pathway']:.2f} | score={c['score']:.2f}"
        )

    print("\n[3] Drug Design & Optimization")
    print("Generated high-potential / low-risk molecules:")
    for g in generated:
        print(
            f"  {g['id']} | efficacy={g['efficacy']:.2f} | tox={g['tox_prob']:.3f} "
            f"| docking={g['docking']:.2f} | pathway={g['pathway']:.2f}"
        )

    plans = stage4_delivery()
    print("\n[4] Smart Drug Delivery")
    for p in plans:
        print(
            f"  {p['patient_id']} | tumor={p['tumor_burden']:.2f} | metastatic={p['metastatic']} "
            f"| immune={p['immune_score']:.2f} | dose={p['recommended_dose_mg']:.1f}mg "
            f"| {p['delivery_mechanism']}"
        )

    write_csv(
        "results_top_driver_genes.csv",
        [{"gene": g, "importance": f"{s:.6f}"} for g, s in stage1.top_driver_genes],
        ["gene", "importance"],
    )
    write_csv(
        "results_top_candidates.csv",
        top_candidates,
        ["id", "mw", "lip", "hbond", "docking", "pathway", "efficacy", "tox_prob", "score"],
    )
    write_csv(
        "results_generated_molecules.csv",
        generated,
        ["id", "efficacy", "tox_prob", "docking", "pathway"],
    )
    write_csv(
        "results_delivery_plan.csv",
        plans,
        ["patient_id", "tumor_burden", "metastatic", "immune_score", "recommended_dose_mg", "delivery_mechanism"],
    )

    print("\nSaved result files:")
    print("- results_top_driver_genes.csv")
    print("- results_top_candidates.csv")
    print("- results_generated_molecules.csv")
    print("- results_delivery_plan.csv")


if __name__ == "__main__":
    main()
