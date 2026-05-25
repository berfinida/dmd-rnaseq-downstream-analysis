#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from sklearn.decomposition import PCA

SAMPLE_COLS = ["DMD1", "DMD2", "WT1", "WT2"]


def load_default_data(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr_path = project_root / "results" / "matrix" / "expression_matrix.tsv"
    summary_path = project_root / "results" / "matrix" / "dmd_vs_wt_summary.tsv"
    expr = pd.read_csv(expr_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")
    return expr, summary


def load_metadata(project_root: Path) -> pd.DataFrame | None:
    metadata_path = project_root / "metadata" / "samples.tsv"
    if not metadata_path.exists():
        return None
    md = pd.read_csv(metadata_path, sep="\t")
    required = {"sample", "group", "condition", "source"}
    if not required.issubset(set(md.columns)):
        return None
    return md


def compute_summary_from_expr(expr: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame()
    out["Name"] = expr["Name"]
    out["DMD_mean_TPM"] = expr[["DMD1", "DMD2"]].mean(axis=1)
    out["WT_mean_TPM"] = expr[["WT1", "WT2"]].mean(axis=1)
    out["log2FC"] = np.log2((out["DMD_mean_TPM"] + 1.0) / (out["WT_mean_TPM"] + 1.0))
    return out


def apply_optional_annotation(summary: pd.DataFrame, project_root: Path) -> tuple[pd.DataFrame, str]:
    mapping_path = project_root / "annotation" / "transcript_to_gene_symbol.tsv"
    if not mapping_path.exists():
        return summary.copy(), "Annotation not available: mapping file not found. Using Ensembl transcript IDs."

    mapping = pd.read_csv(mapping_path, sep="\t")
    required = {"transcript_id", "gene_symbol"}
    if not required.issubset(set(mapping.columns)):
        return summary.copy(), "Annotation not available: mapping columns are invalid. Using Ensembl transcript IDs."

    mapping = mapping[["transcript_id", "gene_symbol"]].drop_duplicates(subset=["transcript_id"])
    merged = summary.merge(mapping, left_on="Name", right_on="transcript_id", how="left").drop(columns=["transcript_id"])
    merged["display_id"] = merged["Name"]
    mask = merged["gene_symbol"].notna() & (merged["gene_symbol"].astype(str).str.strip() != "")
    merged.loc[mask, "display_id"] = merged.loc[mask, "gene_symbol"].astype(str) + " (" + merged.loc[mask, "Name"].astype(str) + ")"
    return merged, "Annotation loaded from annotation/transcript_to_gene_symbol.tsv"


def validate_uploaded_matrix(df: pd.DataFrame) -> tuple[bool, str]:
    required = {"Name", *SAMPLE_COLS}
    missing = required - set(df.columns)
    if missing:
        return False, f"Uploaded matrix is missing required columns: {sorted(missing)}"
    return True, "OK"


def run_pca(expr: pd.DataFrame, metadata: pd.DataFrame | None) -> tuple[pd.DataFrame, np.ndarray]:
    log_expr = np.log2(expr[SAMPLE_COLS] + 1)
    x = log_expr.T.values
    pca = PCA(n_components=2, random_state=42)
    pcs = pca.fit_transform(x)
    explained = pca.explained_variance_ratio_ * 100

    if metadata is not None:
        group_map = dict(zip(metadata["sample"], metadata["group"]))
    else:
        group_map = {"DMD1": "DMD", "DMD2": "DMD", "WT1": "WT", "WT2": "WT"}

    pca_df = pd.DataFrame(
        {
            "sample": SAMPLE_COLS,
            "group": [group_map.get(s, "NA") for s in SAMPLE_COLS],
            "PC1": pcs[:, 0],
            "PC2": pcs[:, 1],
        }
    )
    return pca_df, explained


def main() -> None:
    st.set_page_config(page_title="DMD vs WT Transcriptomics Explorer", page_icon="🧬", layout="wide")
    root = Path(__file__).resolve().parent

    expr_default, summary_default = load_default_data(root)
    metadata = load_metadata(root)

    st.title("DMD vs WT Transcriptomics Explorer")
    st.info("TPM-based descriptive exploratory dashboard. Not formal differential expression analysis.")

    st.sidebar.title("Navigation")
    section = st.sidebar.radio(
        "Go to section",
        [
            "Overview & Data Provenance",
            "Upload Mode",
            "Data Quality & Context",
            "Transcript-Level Explorer",
            "Top Transcript Explorer",
            "Volcano-like Explorer",
            "Heatmap Explorer",
            "PCA (Sample-Level Exploration)",
            "Figure Gallery",
            "Reproducibility",
            "Limitations",
        ],
    )

    st.sidebar.header("Global Filters")
    min_abs_log2fc = st.sidebar.slider("Minimum absolute log2FC", 0.0, 10.0, 2.0, 0.1)
    top_n = st.sidebar.slider("Top N", 5, 200, 25, 5)

    # Upload mode input
    uploaded = st.sidebar.file_uploader("Optional upload: expression_matrix.tsv", type=["tsv", "txt"])
    if uploaded is not None:
        uploaded_expr = pd.read_csv(uploaded, sep="\t")
        ok, msg = validate_uploaded_matrix(uploaded_expr)
        if ok:
            expr = uploaded_expr.copy()
            summary = compute_summary_from_expr(expr)
            st.sidebar.success("Upload validated. Using uploaded matrix for exploratory views.")
        else:
            st.sidebar.error(msg)
            expr = expr_default.copy()
            summary = summary_default.copy()
    else:
        expr = expr_default.copy()
        summary = summary_default.copy()

    # Optional annotation merge
    summary_annot, annotation_msg = apply_optional_annotation(summary, root)
    st.sidebar.caption(annotation_msg)

    if section == "Overview & Data Provenance":
        st.header("Project Overview")
        st.write("Exploratory transcript-level TPM dashboard for DMD vs WT using processed matrix outputs.")

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Number of transcripts", f"{len(expr):,}")
        c2.metric("Number of samples", str(len(SAMPLE_COLS)))
        c3.metric("Max positive log2FC", f"{summary['log2FC'].max():.3f}")
        c4.metric("Max negative log2FC", f"{summary['log2FC'].min():.3f}")

        st.subheader("Data provenance")
        st.markdown(
            "- Source dataset: `GSE156496 / SRP278118`\n"
            "- Original raw data: SRA FASTQ files\n"
            "- Original pipeline: `https://github.com/berfinida/pipeline`\n"
            "- This app uses processed TPM/expression summary files"
        )
        st.warning("All outputs are descriptive exploratory views. No p-values or FDR are computed.")

    elif section == "Upload Mode":
        st.header("Interactive Upload Mode")
        st.write("Upload your own `expression_matrix.tsv` to drive transcript explorer, volcano-like, heatmap, and PCA sections.")
        st.code("Required columns: Name, DMD1, DMD2, WT1, WT2", language="text")
        if uploaded is None:
            st.info("No upload active. Default repository matrix is currently in use.")
        else:
            ok, msg = validate_uploaded_matrix(expr)
            if ok:
                st.success("Uploaded matrix is valid and active.")
                st.dataframe(expr.head(20), use_container_width=True)
            else:
                st.error(msg)

    elif section == "Data Quality & Context":
        st.header("Data Quality & Context")
        missing = int(expr[SAMPLE_COLS].isna().sum().sum())
        tpm_min = float(expr[SAMPLE_COLS].min().min())
        tpm_max = float(expr[SAMPLE_COLS].max().max())
        c1, c2, c3 = st.columns(3)
        c1.metric("Matrix dimensions", f"{expr.shape[0]} x {expr.shape[1]}")
        c2.metric("Missing TPM values", str(missing))
        c3.metric("TPM range", f"{tpm_min:.3f} to {tpm_max:.3f}")

        if metadata is not None:
            st.subheader("Sample metadata")
            st.dataframe(metadata, use_container_width=True)
        else:
            st.warning("Metadata file not available at metadata/samples.tsv")

    elif section == "Transcript-Level Explorer":
        st.header("Transcript-Level Explorer")
        top_ids = (
            summary_annot.assign(abs_log2FC=summary_annot["log2FC"].abs())
            .sort_values("abs_log2FC", ascending=False)
            .head(150)["Name"]
            .astype(str)
            .tolist()
        )

        q = st.text_input("Search by Ensembl transcript ID", placeholder="ENSMUST...")
        choice = st.selectbox("Or choose from high-contrast transcripts", ["(none)"] + top_ids)
        selected = q.strip() if q.strip() else (choice if choice != "(none)" else None)

        if selected:
            expr_hit = expr[expr["Name"].astype(str) == selected]
            summary_hit = summary_annot[summary_annot["Name"].astype(str) == selected]
            if expr_hit.empty:
                st.error("Selected transcript is not present in current matrix.")
            else:
                st.dataframe(expr_hit[["Name", *SAMPLE_COLS]], use_container_width=True)
                cols = [c for c in ["Name", "gene_symbol", "DMD_mean_TPM", "WT_mean_TPM", "log2FC"] if c in summary_hit.columns]
                st.dataframe(summary_hit[cols], use_container_width=True)

                values = expr_hit.iloc[0]
                plot_df = pd.DataFrame({"sample": SAMPLE_COLS, "TPM": [float(values[s]) for s in SAMPLE_COLS]})
                if metadata is not None:
                    gmap = dict(zip(metadata["sample"], metadata["group"]))
                    plot_df["group"] = plot_df["sample"].map(gmap)
                else:
                    plot_df["group"] = plot_df["sample"].str.extract(r"(DMD|WT)")

                fig = px.bar(plot_df, x="sample", y="TPM", color="group", title="Selected transcript TPM")
                st.plotly_chart(fig, use_container_width=True)

    elif section == "Top Transcript Explorer":
        st.header("Top Transcript Explorer")
        sort_mode = st.selectbox(
            "Sort transcripts by",
            [
                "highest log2FC",
                "lowest log2FC",
                "highest absolute log2FC",
                "highest DMD mean TPM",
                "highest WT mean TPM",
            ],
        )

        table = summary_annot.copy()
        table["abs_log2FC"] = table["log2FC"].abs()
        table = table[table["abs_log2FC"] >= min_abs_log2fc]

        if sort_mode == "highest log2FC":
            table = table.sort_values("log2FC", ascending=False)
        elif sort_mode == "lowest log2FC":
            table = table.sort_values("log2FC", ascending=True)
        elif sort_mode == "highest absolute log2FC":
            table = table.sort_values("abs_log2FC", ascending=False)
        elif sort_mode == "highest DMD mean TPM":
            table = table.sort_values("DMD_mean_TPM", ascending=False)
        else:
            table = table.sort_values("WT_mean_TPM", ascending=False)

        top_table = table.head(top_n)
        st.dataframe(top_table, use_container_width=True)

        st.download_button("Download filtered transcript table", top_table.to_csv(sep="\t", index=False), "filtered_transcripts.tsv", "text/tab-separated-values")
        top_up = table.sort_values("log2FC", ascending=False).head(top_n)
        top_down = table.sort_values("log2FC", ascending=True).head(top_n)
        st.download_button("Download top upregulated", top_up.to_csv(sep="\t", index=False), "top_upregulated.tsv", "text/tab-separated-values")
        st.download_button("Download top downregulated", top_down.to_csv(sep="\t", index=False), "top_downregulated.tsv", "text/tab-separated-values")
        st.download_button("Download current descriptive summary", summary_annot.to_csv(sep="\t", index=False), "current_summary.tsv", "text/tab-separated-values")

    elif section == "Volcano-like Explorer":
        st.header("Volcano-like Explorer")
        st.caption("Descriptive TPM-based visualization: log2FC vs abs(log2FC).")
        v = summary_annot.copy()
        v["abs_log2FC"] = v["log2FC"].abs()
        v["highlight"] = v["abs_log2FC"] >= min_abs_log2fc
        hover_cols = [c for c in ["Name", "gene_symbol", "DMD_mean_TPM", "WT_mean_TPM"] if c in v.columns]
        fig = px.scatter(v, x="log2FC", y="abs_log2FC", color="highlight", hover_data=hover_cols)
        fig.add_vline(x=0, line_dash="dash", line_color="gray")
        st.plotly_chart(fig, use_container_width=True)

    elif section == "Heatmap Explorer":
        st.header("Heatmap Explorer")
        st.caption("Top transcripts by absolute log2FC with log2(TPM + 1) values.")
        top = summary_annot.assign(abs_log2FC=summary_annot["log2FC"].abs()).sort_values("abs_log2FC", ascending=False).head(top_n)
        ids = top["Name"].tolist()
        mat = expr[expr["Name"].isin(ids)].set_index("Name").reindex(ids)
        heat = np.log2(mat[SAMPLE_COLS] + 1)
        fig = px.imshow(heat, labels=dict(x="Sample", y="Transcript ID", color="log2(TPM+1)"), aspect="auto")
        st.plotly_chart(fig, use_container_width=True)

    elif section == "PCA (Sample-Level Exploration)":
        st.header("PCA (Sample-Level Exploration)")
        st.markdown("PCA uses log2(TPM + 1) transformed expression values. This is exploratory and not inferential.")
        pca_df, explained = run_pca(expr, metadata)
        fig = px.scatter(pca_df, x="PC1", y="PC2", color="group", text="sample")
        fig.update_layout(xaxis_title=f"PC1 ({explained[0]:.2f}% variance)", yaxis_title=f"PC2 ({explained[1]:.2f}% variance)")
        st.plotly_chart(fig, use_container_width=True)
        st.warning("Small sample size limits interpretation of sample separation.")

    elif section == "Figure Gallery":
        st.header("Figure Gallery")
        fig_dir = root / "downstream_analysis" / "figures"
        items = [
            ("Volcano-like plot", "volcano_like_plot.png", "Descriptive contrast visualization."),
            ("Heatmap", "top_transcripts_heatmap.png", "Top transcripts exploratory heatmap."),
            ("PCA plot", "pca_plot.png", "Sample-level exploratory PCA."),
            ("Summary panel", "summary_panel.png", "Combined descriptive visual summary."),
        ]
        for title, filename, cap in items:
            path = fig_dir / filename
            st.subheader(title)
            st.caption(cap)
            if path.exists():
                st.image(str(path), use_container_width=True)
            else:
                st.info(f"Not available: {path}")

    elif section == "Reproducibility":
        st.header("Reproducibility")
        st.code("pip install -r requirements.txt\nstreamlit run app.py", language="bash")

    elif section == "Limitations":
        st.header("Scientific Limitations")
        st.markdown(
            "- Small sample size (2 DMD, 2 WT)\n"
            "- TPM-based descriptive exploration\n"
            "- No p-values/FDR\n"
            "- Ensembl transcript IDs remain unless real mapping is provided\n"
            "- Not for clinical interpretation"
        )


if __name__ == "__main__":
    main()
