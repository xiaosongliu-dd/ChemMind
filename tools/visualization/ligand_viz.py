from __future__ import annotations

import base64
import html
import json
import tempfile
from pathlib import Path


_VALID_PLOTS = {"grid", "table", "scatter", "umap"}


def run_ligand_viz(
    smiles: list[str],
    properties: dict[str, list] | None = None,
    plots: list[str] | None = None,
    scatter_x: str | None = None,
    scatter_y: str | None = None,
    color_by: str | None = None,
    title: str = "ChemMind ligand report",
    output_dir: str | None = None,
) -> dict:
    """
    Render a self-contained HTML report for a ligand set with structure grid,
    property table, scatter plot, and optional UMAP of Morgan fingerprints.

    The HTML opens in any browser. Designed as a "publish results" step at the
    end of an agentic workflow — chains naturally off ligand_filter, rdkit_enum,
    gnina, boltz2_affinity, etc.
    """
    if plots is None:
        plots = ["grid", "table"]
    unknown = set(plots) - _VALID_PLOTS
    if unknown:
        raise ValueError(f"Unknown plot types: {sorted(unknown)}. "
                         f"Valid: {sorted(_VALID_PLOTS)}")

    properties = properties or {}
    _validate_lengths(smiles, properties)

    out = Path(output_dir or tempfile.mkdtemp(prefix="ligand_viz_"))
    out.mkdir(parents=True, exist_ok=True)

    sections: list[str] = []
    generated: list[str] = []
    skipped:   list[dict] = []

    def _record(plot_name, ok, chunk, reason):
        if ok:
            sections.append(chunk)
            generated.append(plot_name)
        else:
            skipped.append({"plot": plot_name, "reason": reason})

    if "grid" in plots:
        _record("grid", *_render_grid(smiles, properties, out))
    if "table" in plots:
        sections.append(_render_table(smiles, properties))
        generated.append("table")
    if "scatter" in plots:
        _record("scatter", *_render_scatter(
            smiles, properties, scatter_x, scatter_y, color_by, out))
    if "umap" in plots:
        _record("umap", *_render_umap(smiles, properties, color_by, out))

    html_path = out / "report.html"
    html_path.write_text(_compose_html(title, len(smiles), sections))
    return {
        "html_path":          str(html_path),
        "output_dir":         str(out),
        "n_compounds":        len(smiles),
        "plots_generated":    [g for g in generated if isinstance(g, str)],
        "properties_indexed": list(properties.keys()),
        "skipped_plots":      [s for s in skipped if isinstance(s, dict)],
    }


# ── input validation ──────────────────────────────────────────────────────────

def _validate_lengths(smiles: list[str], properties: dict) -> None:
    n = len(smiles)
    for k, v in properties.items():
        if len(v) != n:
            raise ValueError(
                f"Property '{k}' has length {len(v)} but smiles has {n}"
            )


# ── grid (RDKit, optionally mols2grid) ────────────────────────────────────────

def _render_grid(smiles, properties, out_dir: Path) -> tuple[bool, str, str]:
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
    except ImportError as e:
        return False, "", f"rdkit unavailable: {e}"

    try:
        import mols2grid
        ok, chunk, reason = _render_grid_mols2grid(smiles, properties, mols2grid, out_dir)
        if ok:
            return ok, chunk, reason
    except ImportError:
        pass
    except Exception:
        pass  # mols2grid API drift — fall through to static RDKit grid

    mols = [Chem.MolFromSmiles(s) for s in smiles]
    valid = [(s, m) for s, m in zip(smiles, mols) if m is not None]
    if not valid:
        return False, "", "no valid molecules to render"

    n_show = min(len(valid), 100)
    legends = [s if len(s) <= 30 else s[:27] + "..." for s, _ in valid[:n_show]]
    img = Draw.MolsToGridImage(
        [m for _, m in valid[:n_show]],
        molsPerRow=4,
        subImgSize=(250, 200),
        legends=legends,
        useSVG=False,
    )
    png_path = out_dir / "grid.png"
    if hasattr(img, "save"):
        img.save(str(png_path))
    else:
        png_path.write_bytes(img.data)

    b64 = base64.b64encode(png_path.read_bytes()).decode()
    chunk = (
        f'<section><h2>Structure grid '
        f'<small>({n_show} of {len(smiles)})</small></h2>'
        f'<img src="data:image/png;base64,{b64}" '
        f'style="max-width:100%;border:1px solid #ddd"></section>'
    )
    return True, chunk, ""


def _render_grid_mols2grid(smiles, properties, mols2grid, out_dir: Path):
    import pandas as pd

    df = pd.DataFrame({"SMILES": smiles, **properties})
    grid = mols2grid.MolGrid(
        df,
        smiles_col="SMILES",
        prerender=True,
        useSVG=True,
        size=(180, 130),
    )
    grid_path = out_dir / "grid.html"
    grid.save(str(grid_path))
    chunk = (
        '<section><h2>Structure grid '
        '<small>(mols2grid — interactive: filter, sort, select)</small></h2>'
        f'<iframe src="{grid_path.name}" '
        'style="width:100%;height:720px;border:1px solid #ddd"></iframe>'
        '</section>'
    )
    return True, chunk, ""


# ── property table ────────────────────────────────────────────────────────────

def _render_table(smiles, properties) -> str:
    import pandas as pd

    df = pd.DataFrame({"SMILES": smiles, **properties})
    table = df.to_html(index=False, classes="ligand-table",
                       table_id="ligand-table", border=0, escape=True)
    cdn = (
        '<link rel="stylesheet" '
        'href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">'
        '<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>'
        '<script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>'
        '<script>$(function(){$("#ligand-table").DataTable({pageLength:25});});</script>'
    )
    return f'<section><h2>Property table</h2>{cdn}{table}</section>'


# ── scatter plot ──────────────────────────────────────────────────────────────

def _render_scatter(smiles, properties, x_key, y_key, color_by, out_dir
                   ) -> tuple[bool, str, str]:
    if not properties:
        return False, "", "no properties supplied — cannot scatter"
    keys = list(properties.keys())
    x_key = x_key or (keys[0] if len(keys) >= 1 else None)
    y_key = y_key or (keys[1] if len(keys) >= 2 else None)
    if x_key is None or y_key is None or x_key == y_key:
        return False, "", f"need two distinct numeric properties; got x={x_key} y={y_key}"
    if x_key not in properties or y_key not in properties:
        return False, "", f"x/y key not in properties: {x_key}, {y_key}"

    try:
        import plotly.express as px
        import pandas as pd
    except ImportError as e:
        return _render_scatter_mpl(smiles, properties, x_key, y_key, color_by, out_dir)

    df = pd.DataFrame({"SMILES": smiles, **properties})
    df["__row_idx"] = list(range(len(df)))
    color = color_by if color_by in (properties or {}) else None
    fig = px.scatter(df, x=x_key, y=y_key, color=color,
                     custom_data=["__row_idx"])
    fig.update_layout(margin=dict(l=40, r=20, t=40, b=40), height=500)
    fig.update_traces(hoverinfo="none", hovertemplate=None)
    div_id = "ligand-scatter"
    plotly_html = fig.to_html(include_plotlyjs="cdn", full_html=False,
                              div_id=div_id)
    svgs = [_smiles_to_svg(s) for s in smiles]
    overlay = _hover_overlay_js(div_id, svgs)
    chunk = (
        f'<section><h2>{html.escape(y_key)} vs {html.escape(x_key)} '
        f'<small>(hover for structure)</small></h2>'
        f'{plotly_html}{overlay}</section>'
    )
    return True, chunk, ""


def _render_scatter_mpl(smiles, properties, x_key, y_key, color_by, out_dir
                       ) -> tuple[bool, str, str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        return False, "", f"neither plotly nor matplotlib available: {e}"

    fig, ax = plt.subplots(figsize=(7, 5))
    c = properties.get(color_by) if color_by in (properties or {}) else None
    sc = ax.scatter(properties[x_key], properties[y_key], c=c, alpha=0.8)
    if c is not None:
        plt.colorbar(sc, label=color_by)
    ax.set_xlabel(x_key)
    ax.set_ylabel(y_key)
    ax.set_title(f"{y_key} vs {x_key}")
    png = out_dir / "scatter.png"
    fig.savefig(png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    b64 = base64.b64encode(png.read_bytes()).decode()
    chunk = (
        f'<section><h2>{html.escape(y_key)} vs {html.escape(x_key)}</h2>'
        f'<img src="data:image/png;base64,{b64}" '
        f'style="max-width:700px;border:1px solid #ddd"></section>'
    )
    return True, chunk, ""


# ── UMAP of Morgan fingerprints ───────────────────────────────────────────────

def _render_umap(smiles, properties, color_by, out_dir) -> tuple[bool, str, str]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import numpy as np
    except ImportError as e:
        return False, "", f"rdkit/numpy unavailable: {e}"
    try:
        import umap
    except ImportError:
        return False, "", "umap-learn not installed (pip install umap-learn)"

    mols = [Chem.MolFromSmiles(s) for s in smiles]
    fps = []
    keep_idx = []
    for i, m in enumerate(mols):
        if m is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
        fps.append(np.array(fp))
        keep_idx.append(i)
    if len(fps) < 4:
        return False, "", f"need at least 4 valid molecules; got {len(fps)}"

    arr = np.stack(fps)
    n_neighbors = min(15, max(2, len(fps) - 1))
    reducer = umap.UMAP(n_neighbors=n_neighbors, n_components=2,
                        metric="jaccard", random_state=42)
    coords = reducer.fit_transform(arr)

    color_vals = None
    if color_by and color_by in (properties or {}):
        color_vals = [properties[color_by][i] for i in keep_idx]

    try:
        import plotly.express as px
        import pandas as pd
        df = pd.DataFrame({
            "UMAP-1": coords[:, 0],
            "UMAP-2": coords[:, 1],
            "SMILES": [smiles[i] for i in keep_idx],
            "__row_idx": list(range(len(keep_idx))),
        })
        if color_vals is not None:
            df[color_by] = color_vals
        fig = px.scatter(df, x="UMAP-1", y="UMAP-2",
                         color=color_by if color_vals is not None else None,
                         custom_data=["__row_idx"])
        fig.update_layout(height=500)
        fig.update_traces(hoverinfo="none", hovertemplate=None)
        div_id = "ligand-umap"
        plot_html = fig.to_html(include_plotlyjs="cdn", full_html=False,
                                div_id=div_id)
        svgs = [_smiles_to_svg(smiles[i]) for i in keep_idx]
        plot_html += _hover_overlay_js(div_id, svgs)
    except ImportError:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 5))
        sc = ax.scatter(coords[:, 0], coords[:, 1], c=color_vals, alpha=0.8)
        if color_vals is not None:
            plt.colorbar(sc, label=color_by)
        ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
        png = out_dir / "umap.png"
        fig.savefig(png, dpi=120, bbox_inches="tight")
        plt.close(fig)
        b64 = base64.b64encode(png.read_bytes()).decode()
        plot_html = (
            f'<img src="data:image/png;base64,{b64}" '
            f'style="max-width:700px;border:1px solid #ddd">'
        )

    chunk = (
        f'<section><h2>Chemical-space UMAP '
        f'<small>(Morgan FP, Jaccard distance)</small></h2>{plot_html}</section>'
    )
    return True, chunk, ""


# ── hover-over-structure overlay (Plotly) ────────────────────────────────────

def _smiles_to_svg(smi: str, w: int = 220, h: int = 170) -> str:
    """Render a SMILES to a clean inline-embeddable SVG string. Returns "" on
    failure. Strips the leading <?xml?> declaration so the SVG behaves as HTML
    when assigned via innerHTML."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        return ""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return ""
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    if svg.startswith("<?xml"):
        svg = svg.split("?>", 1)[1].lstrip()
    return svg


def _hover_overlay_js(div_id: str, svgs: list[str]) -> str:
    """Return a <script> block that wires plotly_hover on `div_id` to display
    the matching SVG (looked up by customdata[0]) in a floating tooltip."""
    svgs_json = json.dumps(svgs)
    return f"""
<script>
(function() {{
    const SVGS = {svgs_json};
    const div = document.getElementById({json.dumps(div_id)});
    if (!div) return;
    let tip = document.getElementById("ligand-viz-tooltip");
    if (!tip) {{
        tip = document.createElement("div");
        tip.id = "ligand-viz-tooltip";
        tip.style.cssText =
            "position:fixed;pointer-events:none;background:white;"
            + "border:1px solid #888;padding:4px;display:none;z-index:9999;"
            + "box-shadow:0 2px 8px rgba(0,0,0,0.18);border-radius:4px;";
        document.body.appendChild(tip);
    }}
    function show(ev) {{
        if (!ev.points || !ev.points.length) return;
        const cd = ev.points[0].customdata;
        const idx = Array.isArray(cd) ? cd[0] : cd;
        const svg = SVGS[idx];
        if (svg == null) return;
        tip.innerHTML = svg;
        tip.style.display = "block";
        const x = (ev.event && ev.event.clientX) || 0;
        const y = (ev.event && ev.event.clientY) || 0;
        const tw = tip.offsetWidth, th = tip.offsetHeight;
        const vw = window.innerWidth, vh = window.innerHeight;
        tip.style.left = ((x + 14 + tw > vw) ? x - tw - 14 : x + 14) + "px";
        tip.style.top  = ((y + 14 + th > vh) ? y - th - 14 : y + 14) + "px";
    }}
    function hide() {{ tip.style.display = "none"; }}
    function wire() {{
        if (typeof div.on !== "function") {{ setTimeout(wire, 50); return; }}
        div.on("plotly_hover", show);
        div.on("plotly_unhover", hide);
    }}
    wire();
}})();
</script>
"""


# ── HTML composition ──────────────────────────────────────────────────────────

def _compose_html(title: str, n: int, sections: list[str]) -> str:
    style = """
      body { font-family: -apple-system, BlinkMacSystemFont, sans-serif;
             max-width: 1200px; margin: 24px auto; padding: 0 16px; color:#222; }
      h1 { border-bottom: 2px solid #333; padding-bottom: 8px; }
      h2 { margin-top: 32px; color: #444; }
      h2 small { color:#888; font-weight:normal; font-size:0.7em; }
      section { margin-bottom: 32px; }
      table.ligand-table { width:100%; border-collapse: collapse; font-size: 13px; }
      table.ligand-table th, table.ligand-table td {
          padding: 6px 10px; border-bottom: 1px solid #eee; text-align: left; }
      table.ligand-table th { background:#f5f5f5; }
    """
    body = "\n".join(sections) if sections else "<p><em>No plots generated.</em></p>"
    return (
        f'<!DOCTYPE html><html><head><meta charset="utf-8">'
        f'<title>{html.escape(title)}</title>'
        f'<style>{style}</style></head><body>'
        f'<h1>{html.escape(title)}</h1>'
        f'<p><strong>{n}</strong> compounds</p>'
        f'{body}'
        f'</body></html>'
    )
