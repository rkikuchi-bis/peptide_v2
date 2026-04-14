def _normalize_text(x) -> str:
    return str(x).strip().lower()


def calc_target_match_score(target_label: str, title: str) -> float:
    target = _normalize_text(target_label)
    title_l = _normalize_text(title)

    if not target:
        return 0.0

    if target in title_l:
        return 1.0

    target_tokens = [x for x in target.replace("-", " ").split() if x]
    if target_tokens and any(tok in title_l for tok in target_tokens):
        return 0.6

    return 0.0


def calc_ligand_score(title: str) -> float:
    t = _normalize_text(title)
    strong_keywords = ["inhibitor", "ligand", "bound", "complex"]
    if any(k in t for k in strong_keywords):
        return 1.0
    return 0.0


def calc_resolution_score(resolution) -> float:
    try:
        r = float(resolution)
    except Exception:
        return 0.3

    if r <= 2.0:
        return 1.0
    elif r <= 2.5:
        return 0.8
    elif r <= 3.0:
        return 0.6
    elif r <= 3.5:
        return 0.4
    else:
        return 0.2


def calc_method_score(method: str) -> float:
    m = _normalize_text(method)

    if "x-ray" in m or "diffraction" in m:
        return 1.0
    elif "cryo" in m:
        return 0.8
    elif "nmr" in m:
        return 0.5
    elif m.strip():
        return 0.4
    else:
        return 0.3


def calc_domain_score(title: str) -> float:
    t = _normalize_text(title)
    score = 0.5

    if "kinase" in t:
        score += 0.3
    if "inhibitor" in t:
        score += 0.2
    if "complex" in t:
        score += 0.1
    if "extracellular" in t:
        score -= 0.3

    return max(0.0, min(1.0, score))


def calc_structure_priority_score(record: dict, target_label: str) -> dict:
    title = record.get("title", "")
    method = record.get("method", "")
    resolution = record.get("resolution", None)

    target_match_score = calc_target_match_score(target_label, title)
    ligand_score = calc_ligand_score(title)
    resolution_score = calc_resolution_score(resolution)
    method_score = calc_method_score(method)
    domain_score = calc_domain_score(title)

    structure_priority_score = (
        0.35 * target_match_score
        + 0.25 * ligand_score
        + 0.20 * resolution_score
        + 0.10 * method_score
        + 0.10 * domain_score
    )

    out = dict(record)
    out["target_match_score"] = round(target_match_score, 3)
    out["ligand_score"] = round(ligand_score, 3)
    out["resolution_score"] = round(resolution_score, 3)
    out["method_score"] = round(method_score, 3)
    out["domain_score"] = round(domain_score, 3)
    out["structure_priority_score"] = round(structure_priority_score, 3)
    return out
