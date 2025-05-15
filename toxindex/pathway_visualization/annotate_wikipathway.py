import requests
from pathlib import Path
import xmlrpc.client
import matplotlib.colors as mcolors

# -----------------------------------------------------------------------------
# Utility ---------------------------------------------------------------------
# -----------------------------------------------------------------------------

def fetch_gpml(pathway_id: str, dest_dir: str = ".") -> Path:
    """Download the GPML file for *pathway_id* and return the local Path.

    WikiPathways moved away from the legacy *getPathway* web‑service.  The
    current canonical location is the *wikipathways‑assets* bucket.  We keep a
    small list of candidate URLs to stay robust against future moves.
    """
    candidates = [
        # Current production location (May 2025)
        f"https://www.wikipathways.org/wikipathways-assets/pathways/{pathway_id}/{pathway_id}.gpml",
        # Fallback: classic site artefact (rarely needed but harmless)
        f"https://classic.wikipathways.org/downloads/pathway/{pathway_id}.gpml",
    ]

    for url in candidates:
        try:
            r = requests.get(url, timeout = 20)
            if r.status_code == 200 and r.content.strip():
                gpml_path = Path(dest_dir) / f"{pathway_id}.gpml"
                gpml_path.write_bytes(r.content)
                return gpml_path
        except requests.RequestException:
            # Silent continue – try the next candidate URL
            pass

    raise FileNotFoundError(
        f"Unable to retrieve GPML for {pathway_id}. Tried URLs:\n  - " + "\n  - ".join(candidates)
    )

# -----------------------------------------------------------------------------
# Main styling routine ---------------------------------------------------------
# -----------------------------------------------------------------------------

def style_pathway(pathway_id: str, gene_data: dict, output_file: str = "styled_pathway.png") -> None:
    """Render a WikiPathways pathway with custom node sizes and colours.

    Parameters
    ----------
    pathway_id : str
        WikiPathways identifier (e.g. ``'WP3657'``).
    gene_data : dict
        Mapping ``{gene_symbol: {'sum_value': float, 'mean_value': float}}``.
    output_file : str, optional
        Destination PNG file created by PathVisioRPC.
    """

    # 1. Fetch GPML -----------------------------------------------------------
    gpml_path = fetch_gpml(pathway_id)

    # 2. Spawn PathVisio (must be running inside Docker or locally) ----------
    server = xmlrpc.client.ServerProxy("http://localhost:7777")
    server.openPathway(str(gpml_path))

    # 3. Pre‑compute normalisation constants ---------------------------------
    max_mean = max(v["mean_value"] for v in gene_data.values()) or 1.0

    # 4. Iterate over genes ---------------------------------------------------
    for symbol, metrics in gene_data.items():
        sum_val   = metrics["sum_value"]
        mean_val  = metrics["mean_value"]

        # Size: linear scaling (min 60 px, +10 px per unit sum_value)
        width_px  = 60 + int(sum_val * 10)

        # Colour: red channel with alpha proportional to normalised mean_value
        alpha     = mean_val / max_mean
        colour    = mcolors.to_hex(mcolors.to_rgba("red", alpha=alpha))

        try:
            server.setElementColor("GeneProduct", symbol, colour)
            server.setElementProperty("GeneProduct", symbol, "Width", str(width_px))
        except Exception:
            print(f"[warn] Gene symbol '{symbol}' not found in pathway.")

    # 5. Export result --------------------------------------------------------
    server.exportPathwayAsPng(output_file)
    print(f"Styled pathway saved to {output_file}")

# -----------------------------------------------------------------------------
# CLI entry‑point --------------------------------------------------------------
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    example_id = "WP3657"
    example_gene_data = {
        "TP53" : {"sum_value": 5.0, "mean_value": 0.8},
        "BRCA1": {"sum_value": 3.0, "mean_value": 0.6},
        "EGFR" : {"sum_value": 4.0, "mean_value": 0.7},
    }
    style_pathway(example_id, example_gene_data)
