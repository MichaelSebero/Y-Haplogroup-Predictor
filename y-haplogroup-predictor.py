"""
Y-Chromosome Haplogroup Predictor
Aligned to ISOGG 2019 / Underhill et al. 2015 Y-chromosome phylogeny.

SNP Validity Tiers
  Tier 1 — Direct phylogenetically defining SNP; rsID confirmed in dbSNP as the mutation itself.
             Weight: 25 pts.
  Tier 2 — Strong proxy: single homozygous derived state; literature-supported LD with the
             defining SNP; independently replicated in published population studies.
             Weight: 15 pts.
  Tier 3 — Moderate proxy: homozygous derived state accepted; LD support present but weaker
             or less independently replicated. Entries marked NEEDS VERIFICATION should be
             confirmed against a current SNP array manifest or SNP database before deployment.
             Weight: 6 pts.

Y-chromosome genotype rules
  Males carry one Y chromosome; all Y positions are hemizygous.
  Genotyping arrays report Y positions pseudo-homozygously (e.g., A;A, T;T).
  Heterozygous calls (e.g., A;G) on Y positions indicate genotyping error or autosomal
  bleed-through and are NOT accepted as valid matches.

Hierarchical parent-gate scoring
  After raw subclade scoring, any subclade that scores points but whose immediate parent node
  scored zero is treated as an orphan match (likely a false positive) and its contribution to
  the broad haplogroup score is reduced by ORPHAN_GATE_DISCOUNT. This prevents a single deep
  downstream SNP match from generating a spurious broad haplogroup signal.
"""

import os
from collections import defaultdict

# =============================================================================
# SNP Database
# Key structure: haplogroup_snps[broad_hg][subclade][rsid] = [accepted_genotypes]
# All accepted genotypes are homozygous (pseudo-homozygous Y representation).
# Only the derived (mutant) state is listed per SNP.
# =============================================================================
haplogroup_snps = {

    # =========================================================================
    # R1a  (R-M420)
    # Reference: Underhill et al. 2015, ISOGG 2019 tree
    # =========================================================================
    "R1a": {
        # R-M420 — root of all R1a
        "R1a": {
            "rs17222573": ["G;G"],   # Tier 3; M420-region proxy (Underhill 2015 supp.)
            "rs17307677": ["T;T"],   # Tier 3; M420-region proxy
            "rs17306692": ["C;C"],   # Tier 3; M420-region proxy
            "rs17250535": ["T;T"],   # Tier 3; M420-region proxy
        },
        # R-M459 — intermediate node above M198; defines most R1a
        "R1a1": {
            "rs17307105": ["G;G"],   # Tier 3; M459-region proxy
            "rs9786587":  ["T;T"],   # Tier 3; M459-region proxy
            "rs17316227": ["G;G"],   # Tier 3; M459-region proxy
            "rs2534636":  ["A;A"],   # Tier 2; M459-associated (Underhill 2015)
            "rs8179029":  ["G;G"],   # Tier 2; M417/M459-associated (Underhill 2015)
        },
        # R-M198 — R1a1a; the primary R1a clade comprising nearly all R1a lineages
        "R1a1a": {
            "rs2020857":  ["T;T"],   # Tier 2; M198-associated
            "rs1722146":  ["T;T"],   # Tier 2; M198-associated
            "rs17315926": ["C;C"],   # Tier 3; M198-region proxy
            "rs17221601": ["T;T"],   # Tier 3; M198-region proxy
            # rs16981293 REMOVED — this rsID appears in R1b with a different genotype.
            # It cannot be a valid Y-chromosome marker for both haplogroups simultaneously.
        },
        # R-M417 — the critical branching node separating European (Z282) and Asian (Z93) R1a
        "R1a1a1": {
            "rs17316771": ["G;G"],   # Tier 3; M417-region proxy
        },
        # R-M56 — rare, predominantly South Asian subclade
        "R1a1a1a": {
            "rs2032622": ["T;T"],    # Tier 2; M56-associated
        },
        # R-Z282 — all European R1a; parent of Z280, M458, and other European branches
        "R1a1a1b": {
            "rs4988371": ["T;T"],    # Tier 2; Z282-equivalent (Underhill 2015)
            "rs2032659": ["A;A"],    # Tier 2; Z282 branch
        },
        # R-Z280 — Central and Eastern European R1a; largest single European R1a branch
        # Predominant in Poland, Russia, Ukraine, Balkans
        "R1a1a1b1": {
            "rs2032681": ["A;A"],    # Tier 3; Z280-associated — NEEDS VERIFICATION
        },
        # R-M458 — West Slavic R1a; predominant in Poland, Czech Republic, Slovakia
        "R1a1a1b1a": {
            "rs2032645": ["C;C"],    # Tier 3; M458-associated — NEEDS VERIFICATION
        },
        # R-Z93 — Asian R1a; predominant in South Asia, Central Asia, Iran
        "R1a1a1c": {
            "rs2032626": ["G;G"],    # Tier 2; Z93 branch (Underhill 2015)
            "rs2032644": ["C;C"],    # Tier 2; Z93 branch
            "rs2032655": ["G;G"],    # Tier 2; Z93 branch
        },
    },

    # =========================================================================
    # R1b  (R-M343)
    # Reference: Myres et al. 2011, Balaresque et al. 2010, ISOGG 2019 tree
    # =========================================================================
    "R1b": {
        # R-M343 — root of all R1b
        "R1b": {
            "rs9786184": ["A;A"],    # Tier 2; M343-associated (Myres 2011)
            "rs150173":  ["A;A"],    # Tier 2; M343 confirmatory
        },
        # R-M269 — dominant European R1b (~80% of Western European Y chromosomes)
        "R1b1a2": {
            "rs9786882": ["G;G"],    # Tier 3; M269-region proxy
            "rs9786153": ["T;T"],    # Tier 2; M269-associated (Balaresque 2010)
            "rs877756":  ["T;T"],    # Tier 3; M269-region proxy
            "rs2058276": ["G;G"],    # Tier 3; M269-region proxy
            "rs9786140": ["G;G"],    # Tier 3; M269-region proxy
        },
        # R-L11/P310 — backbone ancestor of both P312 and U106; defines most European R1b-M269
        "R1b1a2a": {
            "rs9786048": ["A;A"],    # Tier 3; L11-associated — NEEDS VERIFICATION
        },
        # R-P312/S116 — Western European R1b; parent of L21, DF27, U152
        "R1b1a2a1a1": {
            "rs9786076":  ["T;T"],   # Tier 3; P312-region proxy
            "rs13304168": ["T;T"],   # Tier 3; P312-region proxy
            "rs2082033":  ["T;T"],   # Tier 3; P312-region proxy
            "rs9786283":  ["C;C"],   # Tier 3; P312-region proxy
            "rs9785659":  ["G;G"],   # Tier 3; P312-region proxy
        },
        # R-U106/S21 — Germanic R1b; predominant in Germany, Netherlands, Denmark, England
        "R1b1a2a1a2": {
            "rs34276300": ["A;A"],   # Tier 2; U106-associated
        },
        # R-L21/M529/S145 — Celtic R1b; predominant in Ireland, Scotland, Wales, Brittany
        "R1b1a2a1a1a": {
            "rs11799226": ["A;A"],   # Tier 3; L21-associated — NEEDS VERIFICATION
            # rs16981293 REMOVED — was incorrectly assigned as L21 proxy in the weight system.
            # This rsID conflicts fatally with R1a1a; it has been removed from both haplogroups.
        },
        # R-DF27/S250 — Iberian R1b; predominant in Spain, Portugal, Gascon France
        "R1b1a2a1a1b": {
            "rs17316524": ["T;T"],   # Tier 3; DF27-associated — NEEDS VERIFICATION
        },
        # R-U152/S28 — Italic and Alpine R1b; predominant in Italy, Switzerland, Southern France
        "R1b1a2a1a1c": {
            "rs35277535": ["A;A"],   # Tier 3; U152-associated — NEEDS VERIFICATION
        },
    },

    # =========================================================================
    # I1  (I-M253)  —  Germanic/Nordic/Anglo-Saxon
    # Reference: Rootsi et al. 2004, ISOGG 2019 tree
    # =========================================================================
    "I1": {
        "I1": {
            "rs34046171": ["T;T"],   # TIER 1 — direct M253 defining SNP (dbSNP confirmed)
            "rs17307252": ["G;G"],   # Tier 3; M253-region proxy
            "rs871626":   ["G;G"],   # Tier 3; M253-region proxy
            "rs17221531": ["T;T"],   # Tier 3; M253-region proxy
            "rs9341296":  ["T;T"],   # Tier 2; I1-associated
            "rs13447354": ["A;A"],   # Tier 2; I1-associated
            "rs9341274":  ["G;G"],   # Tier 2; I1-L22 subclade marker
            "rs2032637":  ["G;G"],   # Tier 2; I1-P109 marker
            "rs34626372": ["C;C"],   # Tier 2; I1-DF29 backbone marker (defines majority of I1)
            "rs3912":     ["T;T"],   # Tier 2; I1a marker
        },
    },

    # =========================================================================
    # I2  (I-M438)  —  Balkan / Eastern European
    # Reference: Rootsi et al. 2004, Battaglia et al. 2009, ISOGG 2019 tree
    # =========================================================================
    "I2": {
        # I-M438 — root of I2
        "I2": {
            "rs35547782": ["T;T"],   # Tier 3; M438-region proxy
            "rs17307294": ["G;G"],   # Tier 2; M438-equivalent
        },
        # I-P37.2 — I2a; encompasses the large majority of European I2 lineages
        "I2a": {
            "rs2032674": ["A;A"],    # Tier 3; P37.2-associated — NEEDS VERIFICATION
        },
        # I-M423 — Dinaric / South Slavic subclade (Balkans)
        "I2a1b": {
            "rs2032666": ["T;T"],    # Tier 2; M423-associated (Battaglia 2009)
        },
        # I-M436/P214 — Western European I2 (distinct from I2a)
        "I2b": {
            "rs2032676": ["G;G"],    # Tier 3; I2b-associated — NEEDS VERIFICATION
        },
    },

    # =========================================================================
    # G  (G-M201)  —  Caucasian / Middle Eastern / South European
    # Reference: ISOGG 2019 tree
    # =========================================================================
    "G": {
        # G-M201 — root of G
        "G": {
            "rs2032656": ["T;T"],    # Tier 2; M201-associated
            "rs4988323": ["C;C"],    # Tier 2; G-associated
        },
        # G-P15 — G2; encompasses approximately 90% of all extant G lineages
        "G2": {
            "rs4988324": ["A;A"],    # Tier 3; P15-associated — NEEDS VERIFICATION
        },
        # G-L30/S135 — G2a; dominant subclade in Caucasus, Middle East, Southern Europe
        "G2a": {
            "rs2032698": ["C;C"],    # Tier 3; L30-associated — NEEDS VERIFICATION
        },
    },

    # =========================================================================
    # N  (N-M231)  —  Siberian / Uralic / Finno-Ugric
    # Reference: Rootsi et al. 2007, ISOGG 2019 tree
    # =========================================================================
    "N": {
        # N-M231 — root of N
        "N": {
            "rs4988327": ["G;G"],    # Tier 2; M231-associated (Rootsi 2007)
        },
        # N-M178 — Uralic N; Finnish, Estonian, Saami, Siberian populations
        "N1c": {
            "rs2032688": ["G;G"],    # Tier 2; M178-associated
        },
        # N1c1a1 — dominant Finno-Ugric subclade defined by L550/CTS2929
        "N1c1a1": {
            "rs4988328": ["T;T"],    # Tier 3; L550-associated — NEEDS VERIFICATION
        },
    },
}

# =============================================================================
# Subclade parent map — used for hierarchical gate scoring.
# Key = subclade label, Value = its immediate parent subclade label.
# Any subclade whose parent scored zero has its contribution discounted.
# =============================================================================
SUBCLADE_PARENTS = {
    # R1a hierarchy
    "R1a1":         "R1a",
    "R1a1a":        "R1a1",
    "R1a1a1":       "R1a1a",
    "R1a1a1a":      "R1a1a1",
    "R1a1a1b":      "R1a1a1",
    "R1a1a1b1":     "R1a1a1b",
    "R1a1a1b1a":    "R1a1a1b1",
    "R1a1a1c":      "R1a1a1",
    # R1b hierarchy
    "R1b1a2":       "R1b",
    "R1b1a2a":      "R1b1a2",
    "R1b1a2a1a1":   "R1b1a2a",
    "R1b1a2a1a2":   "R1b1a2a",
    "R1b1a2a1a1a":  "R1b1a2a1a1",
    "R1b1a2a1a1b":  "R1b1a2a1a1",
    "R1b1a2a1a1c":  "R1b1a2a1a1",
    # I2 hierarchy
    "I2a":          "I2",
    "I2a1b":        "I2a",
    "I2b":          "I2",
    # G hierarchy
    "G2":           "G",
    "G2a":          "G2",
    # N hierarchy
    "N1c":          "N",
    "N1c1a1":       "N1c",
}

# Fraction of score credited to broad haplogroup when a subclade is orphaned
# (i.e., it scores but its immediate parent scored zero).
ORPHAN_GATE_DISCOUNT = 0.25

# =============================================================================
# SNP tier weights and tier assignments
# Any SNP not explicitly assigned a tier defaults to Tier 3.
# =============================================================================
TIER_WEIGHTS = {
    1: 25,    # Tier 1 — direct defining SNP
    2: 15,    # Tier 2 — strong proxy; single homozygous derived state; literature-confirmed
    3: 6,     # Tier 3 — moderate proxy; accept with caution; lower discriminatory power
}

SNP_TIERS = {
    # --- Tier 1: direct defining SNPs ---
    "rs34046171": 1,   # I1-M253 (direct dbSNP entry for the M253 mutation)

    # --- Tier 2: strong proxies ---
    # R1a
    "rs2020857":  2,   "rs1722146":  2,   # R1a-M198
    "rs2534636":  2,   "rs8179029":  2,   # R1a-M459/M417
    "rs4988371":  2,   "rs2032659":  2,   # R1a-Z282
    "rs2032626":  2,   "rs2032644":  2,   "rs2032655": 2,  # R1a-Z93
    "rs2032622":  2,                       # R1a-M56
    # R1b
    "rs9786184":  2,   "rs150173":   2,   # R1b-M343
    "rs9786153":  2,                       # R1b-M269
    "rs34276300": 2,                       # R1b-U106
    # I1
    "rs9341296":  2,   "rs13447354": 2,   # I1 general
    "rs9341274":  2,   "rs2032637":  2,   # I1-L22, I1-P109
    "rs34626372": 2,   "rs3912":     2,   # I1-DF29, I1a
    # I2
    "rs17307294": 2,   "rs2032666":  2,   # I2-M438, I2a1b-M423
    # G
    "rs2032656":  2,   "rs4988323":  2,   # G-M201
    # N
    "rs4988327":  2,   "rs2032688":  2,   # N-M231, N1c-M178
}


def get_snp_weight(snp):
    """Return the weight for a SNP based on its tier assignment."""
    tier = SNP_TIERS.get(snp, 3)   # unassigned SNPs are Tier 3
    return TIER_WEIGHTS[tier]


def get_snp_tier(snp):
    """Return the tier number for a SNP (1, 2, or 3)."""
    return SNP_TIERS.get(snp, 3)


# =============================================================================
# File parsing
# =============================================================================

def parse_dna_file(file_path):
    """
    Parse a DNA file in 23andMe or AncestryDNA format.

    Supported formats:
      4-column:  rsid  chr  pos  genotype       (e.g., 23andMe v3/v4/v5)
      5-column:  rsid  chr  pos  allele1  allele2  (e.g., AncestryDNA)

    Heterozygous calls on Y-chromosome positions indicate a genotyping artifact.
    They are parsed and stored, but will not match any accepted (homozygous-only)
    genotype in haplogroup_snps and will therefore never contribute to scores.
    """
    detected_snps = {}
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#") or "rsid" in line.lower() or line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            rsid = parts[0]
            if len(parts) == 4:
                genotype = parts[3]
                if len(genotype) == 2:
                    allele1, allele2 = genotype[0], genotype[1]
                else:
                    continue
            else:
                allele1, allele2 = parts[3], parts[4]

            if allele1 == allele2:
                genotype_formatted = f"{allele1};{allele2}"
            else:
                # Normalise to alphabetical order for consistent matching
                a, b = sorted([allele1, allele2])
                genotype_formatted = f"{a};{b}"

            detected_snps[rsid] = genotype_formatted

    return detected_snps


# =============================================================================
# Scoring
# =============================================================================

def calculate_haplogroups(detected_snps):
    """
    Two-pass scoring:

    Pass 1 — Raw subclade scoring
      For each SNP in each subclade, if the detected genotype exactly matches
      one of the accepted (homozygous derived) states, add the SNP's tier weight
      to that subclade's raw score.
      An empty genotype list is treated as undefined and scores nothing.

    Pass 2 — Hierarchical gate and broad accumulation
      For each subclade with a non-zero raw score, check its parent in
      SUBCLADE_PARENTS. If the parent scored zero, apply ORPHAN_GATE_DISCOUNT
      to the subclade's contribution to the broad haplogroup total. This
      reduces false positives from isolated deep-subclade matches.
    """
    subclade_raw_scores  = defaultdict(int)
    snps_detected_in_haplogroups = defaultdict(list)

    # Pass 1: score each subclade independently
    for broad_hg, subclades in haplogroup_snps.items():
        for subclade, snps_list in subclades.items():
            for snp, genotypes in snps_list.items():
                if snp not in detected_snps:
                    continue
                # BUG FIX: empty genotype list must never score
                if not genotypes:
                    continue
                if detected_snps[snp] in genotypes:
                    subclade_raw_scores[subclade] += get_snp_weight(snp)
                    snps_detected_in_haplogroups[subclade].append(snp)

    # Pass 2: hierarchical gate + accumulate broad scores
    broad_haplogroup_scores = defaultdict(float)
    for broad_hg, subclades in haplogroup_snps.items():
        for subclade in subclades:
            raw = subclade_raw_scores[subclade]
            if raw == 0:
                continue
            parent = SUBCLADE_PARENTS.get(subclade)
            if parent is not None and subclade_raw_scores[parent] == 0:
                # Orphaned match: discount contribution
                contribution = raw * ORPHAN_GATE_DISCOUNT
            else:
                contribution = float(raw)
            broad_haplogroup_scores[broad_hg] += contribution

    # Build specific scores dict (raw, unaffected by gate)
    specific_haplogroup_scores = {}
    for broad_hg, subclades in haplogroup_snps.items():
        for subclade in subclades:
            score = subclade_raw_scores[subclade]
            specific_haplogroup_scores[subclade] = score if score > 0 else "not calculated"

    return broad_haplogroup_scores, specific_haplogroup_scores, snps_detected_in_haplogroups


# =============================================================================
# Prediction
# =============================================================================

def predict_haplogroup(broad_haplogroup_scores):
    """
    Predict the most likely broad haplogroup.
    Returns 'Unknown (low confidence)' if the top score is below threshold.
    Reports close competitors (within 80% of top score) as alternatives.
    """
    if not broad_haplogroup_scores:
        return "Unknown"

    predicted = max(broad_haplogroup_scores, key=broad_haplogroup_scores.get)
    max_score = broad_haplogroup_scores[predicted]

    if max_score < 15:
        return "Unknown (low confidence)"

    close_competitors = [
        hg for hg, score in broad_haplogroup_scores.items()
        if score >= max_score * 0.8 and hg != predicted
    ]

    if close_competitors:
        return f"{predicted} (possible: {', '.join(close_competitors)})"

    return predicted


def get_most_specific_subclade(specific_haplogroup_scores, predicted_broad_haplogroup):
    """
    Among all subclades belonging to the predicted broad haplogroup,
    return the one that is most specific (longest label) with the highest raw score.
    """
    if predicted_broad_haplogroup in ("Unknown", "Unknown (low confidence)"):
        return "Unknown"

    base = predicted_broad_haplogroup.split()[0]
    candidates = {
        sc: score
        for sc, score in specific_haplogroup_scores.items()
        if sc.startswith(base) and isinstance(score, int) and score > 0
    }

    if not candidates:
        return base

    best = max(candidates.items(), key=lambda x: (x[1], len(x[0])))
    return best[0]


# =============================================================================
# Results display
# =============================================================================

BOLD  = "\033[1m"
RESET = "\033[0m"
DIM   = "\033[2m"
YELLOW = "\033[33m"
GREEN  = "\033[32m"
RED    = "\033[31m"
CYAN   = "\033[36m"

TIER_LABELS = {1: f"{GREEN}[T1 Direct]{RESET}", 2: f"{CYAN}[T2 Proxy]{RESET}", 3: f"{DIM}[T3 Proxy]{RESET}"}

# SNPs that are new proxies requiring independent validation
NEEDS_VERIFICATION = {
    "rs2032681", "rs2032645",                            # R1a-Z280, M458
    "rs9786048",                                          # R1b-L11
    "rs11799226", "rs17316524", "rs35277535",            # R1b-L21, DF27, U152
    "rs2032674", "rs2032676",                            # I2a, I2b
    "rs4988324", "rs2032698",                            # G2, G2a
    "rs4988328",                                          # N1c1a1
}


def print_results(predicted_haplogroup, broad_haplogroup_scores,
                  specific_haplogroup_scores, snps_detected_in_haplogroups,
                  detected_snps):
    """Print formatted prediction results with tier annotations and validation flags."""

    print(f"\n{BOLD}Predicted Y Haplogroup:{RESET} {predicted_haplogroup}")

    most_specific = get_most_specific_subclade(specific_haplogroup_scores, predicted_haplogroup)
    if most_specific not in (predicted_haplogroup.split()[0], "Unknown"):
        print(f"{BOLD}Most Specific Subclade: {RESET}{most_specific}")

    # Confidence
    max_score = max(broad_haplogroup_scores.values(), default=0)
    if max_score >= 40:
        conf_str = f"{GREEN}High{RESET}"
    elif max_score >= 20:
        conf_str = f"{YELLOW}Medium{RESET}"
    else:
        conf_str = f"{RED}Low{RESET}"
    print(f"{BOLD}Confidence Level:{RESET} {conf_str}  (top haplogroup score: {max_score:.1f} pts)")

    # Broad scores
    print(f"\n{BOLD}Broad Haplogroup Scores:{RESET}")
    for hg, score in sorted(broad_haplogroup_scores.items(), key=lambda x: x[1], reverse=True):
        print(f"  {hg}: {score:.1f} pts")

    # Subclade breakdown
    print(f"\n{BOLD}Subclade Scores:{RESET}")
    for broad_hg in sorted(broad_haplogroup_scores):
        relevant = [
            (sc, score) for sc, score in specific_haplogroup_scores.items()
            if sc.startswith(broad_hg) and isinstance(score, int) and score > 0
        ]
        if not relevant:
            continue
        print(f"  {broad_hg}:")
        for subclade, score in sorted(relevant, key=lambda x: x[1], reverse=True):
            # Check orphan status
            parent = SUBCLADE_PARENTS.get(subclade)
            orphan_flag = ""
            if parent:
                parent_score = specific_haplogroup_scores.get(parent, "not calculated")
                if parent_score == "not calculated" or parent_score == 0:
                    orphan_flag = f"  {YELLOW}[gate discounted]{RESET}"
            print(f"    {subclade}: {score} pts (raw){orphan_flag}")

    # Detected SNPs with tier annotations
    print(f"\n{BOLD}Detected SNPs by Subclade:{RESET}")
    any_needs_verification = False
    for subclade, snps in snps_detected_in_haplogroups.items():
        if not snps:
            continue
        snp_strings = []
        for snp in snps:
            tier = get_snp_tier(snp)
            label = TIER_LABELS[tier]
            genotype = detected_snps.get(snp, "?")
            needs_v = f" {YELLOW}*{RESET}" if snp in NEEDS_VERIFICATION else ""
            snp_strings.append(f"{snp} ({genotype}) {label}{needs_v}")
            if snp in NEEDS_VERIFICATION:
                any_needs_verification = True
        print(f"  {subclade}:")
        for s in snp_strings:
            print(f"    {s}")

    if any_needs_verification:
        print(f"\n  {YELLOW}*{RESET} Proxy SNP marked NEEDS VERIFICATION — confirm rsID against")
        print( "    current array manifest or SNP database before treating as definitive.")

    # SNP tier legend
    print(f"\n{DIM}Tier legend: [T1] direct defining SNP (25 pts)  "
          f"[T2] strong proxy (15 pts)  [T3] moderate proxy (6 pts){RESET}")


# =============================================================================
# Main
# =============================================================================

def main():
    print(f"{BOLD}Y-Chromosome Haplogroup Predictor{RESET}")
    print("=" * 40)

    dna_file_path = input("Enter the path to your DNA file: ").strip()

    if not os.path.exists(dna_file_path):
        print("Error: The specified file does not exist.")
        return

    try:
        print("\nParsing DNA file...")
        detected_snps = parse_dna_file(dna_file_path)

        if not detected_snps:
            print("Error: No valid SNP data found in the file.")
            return

        print(f"Found {len(detected_snps):,} SNPs.")

        print("Scoring haplogroups...")
        broad_scores, specific_scores, snps_by_subclade = calculate_haplogroups(detected_snps)

        predicted = predict_haplogroup(broad_scores)

        print("\n" + "=" * 60)
        print_results(predicted, broad_scores, specific_scores, snps_by_subclade, detected_snps)

    except Exception as e:
        print(f"\nError processing file: {e}")
        print("Please verify the file is in 23andMe or AncestryDNA tab-separated format.")


if __name__ == "__main__":
    main()
