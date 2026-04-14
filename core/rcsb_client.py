import json
from urllib import request, error

from core.structure_scorer import calc_structure_priority_score


class InMemoryUploadedStructure:
    def __init__(self, data: bytes, name: str):
        self._data = data
        self.name = name

    def getvalue(self):
        return self._data


def http_post_json(url: str, payload: dict, timeout: int = 20) -> dict:
    body = json.dumps(payload).encode("utf-8")
    req = request.Request(
        url,
        data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def http_get_json(url: str, timeout: int = 20) -> dict:
    req = request.Request(url, headers={"Accept": "application/json"})
    with request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def http_get_bytes(url: str, timeout: int = 30) -> bytes:
    req = request.Request(url, headers={"Accept": "*/*"})
    with request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _safe_nested_get(d: dict, path: list, default=None):
    cur = d
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def fetch_rcsb_entry_metadata(pdb_id: str) -> dict:
    pdb_id = str(pdb_id).lower()
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        data = http_get_json(url)
    except Exception:
        return {
            "pdb_id": pdb_id.upper(),
            "title": "",
            "method": "",
            "resolution": None,
            "deposit_date": "",
        }

    title = _safe_nested_get(data, ["struct", "title"], "")
    exptl = data.get("exptl", [])
    method = ""
    if isinstance(exptl, list) and len(exptl) > 0:
        method = exptl[0].get("method", "")

    resolution = None
    rcsb_entry_info = data.get("rcsb_entry_info", {})
    if isinstance(rcsb_entry_info, dict):
        res_list = rcsb_entry_info.get("resolution_combined")
        if isinstance(res_list, list) and len(res_list) > 0:
            resolution = res_list[0]

    deposit_date = ""
    accession_info = data.get("rcsb_accession_info", {})
    if isinstance(accession_info, dict):
        deposit_date = accession_info.get("deposit_date", "")

    return {
        "pdb_id": pdb_id.upper(),
        "title": title,
        "method": method,
        "resolution": resolution,
        "deposit_date": deposit_date,
    }


def search_rcsb_structures(query_text: str, target_label: str, rows: int = 10) -> list[dict]:
    """
    Search RCSB by free text, enrich with entry metadata,
    then add structure priority score and sort descending.
    Handles both compact(string) and object-style result_set items safely.
    """
    q = str(query_text).strip()
    if not q:
        return []

    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": q,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": rows,
            },
            "results_verbosity": "minimal",
            "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc",
                }
            ],
        },
    }

    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    data = http_post_json(search_url, payload)

    result_set = data.get("result_set", [])
    pdb_ids = []

    for item in result_set:
        identifier = None

        if isinstance(item, str):
            identifier = item
        elif isinstance(item, dict):
            identifier = item.get("identifier")

        if identifier:
            pdb_ids.append(str(identifier).upper())

    out = []
    for pdb_id in pdb_ids:
        meta = fetch_rcsb_entry_metadata(pdb_id)
        meta = calc_structure_priority_score(meta, target_label=target_label)
        out.append(meta)

    out = sorted(
        out,
        key=lambda x: x.get("structure_priority_score", 0),
        reverse=True,
    )
    return out


def download_rcsb_mmcif(pdb_id: str) -> tuple[bytes, str]:
    pid = str(pdb_id).lower()
    url = f"https://files.rcsb.org/download/{pid}.cif"
    file_bytes = http_get_bytes(url)
    filename = f"{pid}.cif"
    return file_bytes, filename


def format_resolution(x):
    if x is None or x == "":
        return "-"
    try:
        return f"{float(x):.2f} Å"
    except Exception:
        return str(x)


def build_rcsb_label(record: dict) -> str:
    pdb_id = record.get("pdb_id", "-")
    title = record.get("title", "") or "-"
    resolution = format_resolution(record.get("resolution"))
    score = record.get("structure_priority_score", "-")
    return f"{pdb_id} | score={score} | {title} | {resolution}"
