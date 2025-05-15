import pywikipathways as pw
import xmlrpc.client
import matplotlib.colors as mcolors

def style_pathway(pathway_id, gene_data, output_file="styled_pathway.png"):
    """
    pathway_id: str (e.g., 'WP3657')
    gene_data: dict of form { "GENE_SYMBOL": {"sum_value": float, "mean_value": float}, ... }
    """

    # Connect to PathVisio RPC (PathVisio should be running)
    server = xmlrpc.client.ServerProxy("http://localhost:7777")

    # Download the pathway GPML file
    gpml_url = pw.getPathwayInfo(pathway_id)['gpml']
    gpml_filename = f"{pathway_id}.gpml"
    with open(gpml_filename, "wb") as f:
        f.write(pw.core._download(gpml_url))  # internal helper

    # Open pathway in PathVisio
    server.openPathway(gpml_filename)

    # For each gene, apply box size and color
    for gene_symbol, values in gene_data.items():
        sum_value = values["sum_value"]
        mean_value = values["mean_value"]

        # Normalize mean_value to [0,1] for shading
        norm_color = mcolors.to_hex(mcolors.to_rgba("red", alpha=mean_value / max(1e-6, max(v["mean_value"] for v in gene_data.values()))))

        # Use sum_value to scale width (arbitrary mapping)
        box_width = 60 + int(sum_value * 10)  # example mapping

        # Apply style
        try:
            server.setElementColor("GeneProduct", gene_symbol, norm_color)
            server.setElementProperty("GeneProduct", gene_symbol, "Width", str(box_width))
        except Exception:
            print(f"Warning: Gene {gene_symbol} not found in pathway")

    # Export the styled image
    server.exportPathwayAsPng(output_file)
    print(f"Styled pathway saved to {output_file}")
