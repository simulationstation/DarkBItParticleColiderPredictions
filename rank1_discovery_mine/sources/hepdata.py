"""
HEPData source module.

HEPData (https://www.hepdata.net) is the primary source for published
high-energy physics data tables.

API Documentation: https://www.hepdata.net/api
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
from urllib.parse import quote_plus

logger = logging.getLogger(__name__)

HEPDATA_BASE = "https://www.hepdata.net"
HEPDATA_SEARCH_URL = "https://www.hepdata.net/search"


def plan_queries(candidate: "Candidate") -> List[str]:
    """
    Plan HEPData search queries for a candidate.

    Returns list of query URLs that would be executed.
    """
    queries = []

    # Query by search terms
    for term in candidate.search_terms:
        encoded = quote_plus(term)
        queries.append(f"{HEPDATA_SEARCH_URL}?q={encoded}&format=json")

    # Query by known HEPData IDs (INSPIRE IDs)
    for hep_id in candidate.hepdata_ids:
        queries.append(f"{HEPDATA_BASE}/record/ins{hep_id}?format=json")

    # Query by collaboration + states
    if candidate.collaboration:
        for state in candidate.states:
            term = f"{candidate.collaboration} {state}"
            encoded = quote_plus(term)
            queries.append(f"{HEPDATA_SEARCH_URL}?q={encoded}&format=json")

    return queries


def execute_search(candidate: "Candidate") -> List[Dict[str, Any]]:
    """
    Execute HEPData searches and return discovered URLs.

    NOTE: This actually makes network requests. Only call when --execute is set.

    Returns list of dicts with keys: url, type, record_id, table_id, description

    HEPData ID formats:
    - HEPData internal ID: numeric (e.g., 89271) -> /record/89271
    - INSPIRE ID: with 'ins' prefix (e.g., ins1728691) -> /record/ins1728691

    Search results return HEPData internal IDs (no 'ins' prefix).
    """
    import requests

    discovered = []
    queries = plan_queries(candidate)

    for query_url in queries:
        try:
            logger.info(f"HEPData query: {query_url}")
            response = requests.get(query_url, timeout=30)
            response.raise_for_status()

            data = response.json()

            # Handle search results
            if "results" in data:
                for result in data["results"]:
                    record_id = result.get("recid")
                    if record_id:
                        # Search returns HEPData internal ID, NOT INSPIRE ID
                        # Use record_id directly without 'ins' prefix
                        discovered.append({
                            "url": f"{HEPDATA_BASE}/record/{record_id}?format=json",
                            "type": "hepdata_record",
                            "record_id": record_id,
                            "title": result.get("title", ""),
                            "collaboration": result.get("collaboration", ""),
                        })

            # Handle direct record lookup (from INSPIRE ID query)
            elif "record" in data:
                record = data["record"]
                record_id = record.get("recid")
                tables = record.get("tables", [])
                for table in tables:
                    table_id = table.get("id")
                    discovered.append({
                        "url": f"{HEPDATA_BASE}/download/table/{record_id}/{table_id}?format=csv",
                        "type": "hepdata_table",
                        "record_id": record_id,
                        "table_id": table_id,
                        "description": table.get("description", ""),
                    })

        except Exception as e:
            logger.warning(f"HEPData query failed: {e}")

    # Deduplicate by URL
    seen = set()
    unique = []
    for item in discovered:
        if item["url"] not in seen:
            seen.add(item["url"])
            unique.append(item)

    return unique


def download(url: str, dest_dir: Path, timeout: int = 60) -> List[str]:
    """
    Download HEPData table(s) to destination directory.

    Returns list of downloaded file paths.

    URL formats supported:
    - /record/{id}?format=json - download full record JSON + all tables as CSV
    - /download/table/{inspire_id}/{table_name}/{version}/csv - download single table

    The record_id can be either HEPData internal ID or INSPIRE ID (with ins prefix).
    """
    import requests
    from urllib.parse import quote

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    try:
        # Determine if this is a table download or record URL
        if "/download/table/" in url:
            # Direct table CSV download
            logger.debug(f"Downloading table CSV: {url}")
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()

            # Extract filename from URL
            path_part = url.split("?")[0].rstrip("/csv")
            parts = path_part.split("/")
            # Handle /download/table/{inspire_id}/{table_name}/{version}
            table_name = parts[-2] if len(parts) > 2 else parts[-1]
            inspire_id = parts[-3] if len(parts) > 3 else "unknown"
            # Clean table name for filename
            safe_name = table_name.replace(" ", "_").replace("/", "_")
            filename = f"hepdata_{inspire_id}_{safe_name}.csv"
            filepath = dest_dir / filename

            with open(filepath, 'wb') as f:
                f.write(response.content)
            downloaded.append(str(filepath))
            logger.info(f"Downloaded: {filename}")

        elif "/record/" in url:
            # Full record - fetch JSON first
            # Ensure we request JSON format
            if "format=json" not in url:
                url = url + ("&" if "?" in url else "?") + "format=json"

            logger.debug(f"Fetching record: {url}")
            response = requests.get(url, timeout=timeout)
            response.raise_for_status()

            # Extract record_id from URL: /record/{id}?format=json
            path_part = url.split("?")[0]
            record_id = path_part.split("/")[-1]
            filename = f"hepdata_{record_id}.json"
            filepath = dest_dir / filename

            with open(filepath, 'wb') as f:
                f.write(response.content)
            downloaded.append(str(filepath))
            logger.info(f"Downloaded: {filename}")

            # Get the INSPIRE ID and version from the response
            data = response.json()
            inspire_id = data.get("record", {}).get("inspire_id")
            version = data.get("version", 1)

            # Download individual tables as CSV
            # Tables are in data_tables at top level, not in record
            tables = data.get("data_tables", [])
            if tables:
                logger.info(f"Record has {len(tables)} tables")
                for table in tables:
                    table_name = table.get("name")
                    if table_name and inspire_id:
                        try:
                            # URL format: /download/table/{inspire_id}/{table_name}/{version}/csv
                            encoded_name = quote(table_name, safe='')
                            table_url = f"{HEPDATA_BASE}/download/table/ins{inspire_id}/{encoded_name}/{version}/csv"
                            table_files = download(table_url, dest_dir, timeout)
                            downloaded.extend(table_files)
                        except Exception as te:
                            logger.warning(f"Failed to download table '{table_name}': {te}")

    except requests.exceptions.HTTPError as e:
        logger.error(f"HEPData HTTP error ({e.response.status_code}): {url}")
    except Exception as e:
        logger.error(f"HEPData download failed: {e}")

    return downloaded


def get_record_metadata(record_id: str) -> Optional[Dict]:
    """
    Get metadata for a HEPData record.

    Args:
        record_id: HEPData record ID (numeric) or INSPIRE ID (with 'ins' prefix)

    NOTE: Makes network request.
    """
    import requests

    try:
        # Handle both HEPData internal ID and INSPIRE ID
        if str(record_id).startswith("ins"):
            url = f"{HEPDATA_BASE}/record/{record_id}?format=json"
        else:
            url = f"{HEPDATA_BASE}/record/{record_id}?format=json"
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        logger.error(f"Failed to get HEPData metadata: {e}")
        return None
