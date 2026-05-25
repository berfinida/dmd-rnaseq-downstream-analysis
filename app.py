#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from sklearn.decomposition import PCA

SAMPLE_COLS = ["DMD1", "DMD2", "WT1", "WT2"]
GROUP_MAP = {"DMD1": "DMD", "DMD2": "DMD", "WT1": "WT", "WT2": "WT"}


def load_data(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr_path = project_root / "results" / "matrix" / "expression_matrix.tsv"
    summary_path = project_root / "results" / "matrix" / "dmd_vs_wt_summary.tsv"
    expr = pd.read_csv(expr_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")
    return expr, summary


def run_pca(expr: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    log_expr = np.log2(expr[SAMPLE_COLS] + 1)
    x = log_expr.T.values
    pca = PCA(n_components=2, random_state=42)
    pcs = pca.fit_transform(x)
    explained = pca.explained_variance_ratio_ * 100
    pca_df = pd.DataFrame(
        {
            "sample": SAMPLE_COLS,
            "group": [GROUP_MAP[s] for s in SAMPLE_COLS],
            "PC1": pcs[:, 0],
            "PC2": pcs[:, 1],
        }
    )
    return pca_df, explained


def main() -> None:
    st.set_page_config(page_title="DMD vs WT Transcriptomics Explorer", page_icon="🧬", layout="wide")

    project_root = Path(__file__).resolve().parent
    expr, summary = load_data(project_root)

    required_expr = {"Name", *SAMPLE_COLS}
    required_summary = {"Name", "DMD_mean_TPM", "WT_mean_TPM", "log2FC"}
    if not required_expr.issubset(set(expr.columns)):
        st.error("Expression matrix is missing required columns.")
        return
    if not required_summary.issubset(set(summary.columns)):
        st.error("Summary table is missing required columns.")
        return

    st.title("DMD vs WT Transcriptomics Explorer")
    st.info(
        "This dashboard is for TPM-based descriptive exploratory analysis. "
        "It does not perform formal statistical differential expression testing."
    )

    st.sidebar.title("Navigation")
    section = st.sidebar.radio(
        "Go to section",
        [
            "Overview & Data Provenance",
            "Data Quality & Context",
            "Transcript-Level Explorer",
            "Top Transcript Explorer",
            "PCA (Sample-Level Exploration)",
            "Figure Gallery",
            "Reproducibility",
            "Limitations",
        ],
    )

    st.sidebar.header("Top-Table Filters")
    min_abs_log2fc = st.sidebar.slider("Minimum absolute log2FC", 0.0, 10.0, 2.0, 0.1)
    top_n = st.sidebar.slider("Top N", 5, 200, 25, 5)

    if section == "Overview & Data Provenance":
        st.header("Project Overview")
        st.write(
            "This app explores transcript-level TPM summaries for DMD vs WT samples "
            "using processed matrix outputs from an upstream RNA-seq pipeline."
        )

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Number of transcripts", f"{len(expr):,}")
        c2.metric("Number of samples", str(len(SAMPLE_COLS)))
        c3.metric("Max positive log2FC", f"{summary['log2FC'].max():.3f}")
        c4.metric("Max negative log2FC", f"{summary['log2FC'].min():.3f}")

        st.subheader("Data Provenance")
        st.markdown(
            "- Source dataset: `GSE156496` / `SRP278118`\n"
            "- Original raw data: SRA FASTQ files\n"
            "- Original pipeline: Nextflow RNA-seq pipeline (`https://github.com/berfinida/pipeline`)\n"
            "- This app uses processed files: `results/matrix/expression_matrix.tsv` and `results/matrix/dmd_vs_wt_summary.tsv`"
        )

        st.warning(
            "Gene-symbol annotation is not applied in this app by default. "
            "Ensembl transcript IDs are shown unless an external transcript-to-symbol mapping is provided."
        )

    elif section == "Data Quality & Context":
        st.header("Data Quality & Context")

        missing_expr = int(expr[SAMPLE_COLS].isna().sum().sum())
        tpm_min = float(expr[SAMPLE_COLS].min().min())
        tpm_max = float(expr[SAMPLE_COLS].max().max())

        c1, c2, c3 = st.columns(3)
        c1.metric("Matrix dimensions", f"{expr.shape[0]} x {expr.shape[1]}")
        c2.metric("Missing TPM values", str(missing_expr))
        c3.metric("TPM range", f"{tpm_min:.3f} to {tpm_max:.3f}")

        st.markdown(
            "**Sample names:** `DMD1`, `DMD2`, `WT1`, `WT2`  \n"
            "**Group labels:** `DMD` = (`DMD1`, `DMD2`), `WT` = (`WT1`, `WT2`)"
        )

        st.subheader("Expression Matrix Preview")
        st.dataframe(expr.head(30), use_container_width=True)

    elif section == "Transcript-Level Explorer":
        st.header("Transcript-Level Explorer")
        st.caption("Search using Ensembl transcript IDs present in the matrix.")

        ranked_ids = (
            summary.assign(abs_log2FC=summary["log2FC"].abs())
            .sort_values("abs_log2FC", ascending=False)
            .head(150)["Name"]
            .astype(str)
            .tolist()
        )

        query = st.text_input("Search by Ensembl transcript ID", placeholder="ENSMUST...")
        drop_choice = st.selectbox("Or choose from high-contrast transcripts", ["(none)"] + ranked_ids)
        selected = query.strip() if query.strip() else (drop_choice if drop_choice != "(none)" else None)

        if selected:
            expr_hit = expr[expr["Name"].astype(str) == selected]
            summary_hit = summary[summary["Name"].astype(str) == selected]

            if expr_hit.empty:
                st.error("Selected transcript ID was not found in the expression matrix.")
            else:
                st.subheader(f"Transcript: {selected}")
                st.dataframe(expr_hit[["Name", *SAMPLE_COLS]], use_container_width=True)

                if summary_hit.empty:
                    st.warning("Transcript was found in expression matrix but not in summary table.")
                else:
                    st.dataframe(
                        summary_hit[["Name", "DMD_mean_TPM", "WT_mean_TPM", "log2FC"]],
                        use_container_width=True,
                    )

                vals = expr_hit.iloc[0]
                plot_df = pd.DataFrame(
                    {
                        "sample": SAMPLE_COLS,
                        "TPM": [float(vals[s]) for s in SAMPLE_COLS],
                        "group": [GROUP_MAP[s] for s in SAMPLE_COLS],
                    }
                )
                fig = px.bar(
                    plot_df,
                    x="sample",
                    y="TPM",
                    color="group",
                    title="Selected transcript TPM across samples",
                    color_discrete_map={"DMD": "#D1495B", "WT": "#2E86AB"},
                )
                st.plotly_chart(fig, use_container_width=True)

    elif section == "Top Transcript Explorer":
        st.header("Top Transcript Explorer")
        st.caption("Sort and filter transcript-level summary values for exploratory inspection.")

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

        table = summary.copy()
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
        st.download_button(
            "Download filtered table (TSV)",
            data=top_table.to_csv(sep="\t", index=False),
            file_name="top_transcript_explorer.tsv",
            mime="text/tab-separated-values",
        )

    elif section == "PCA (Sample-Level Exploration)":
        st.header("PCA (Sample-Level Exploratory Visualization)")
        st.markdown(
            "PCA is computed on `log2(TPM + 1)` transformed transcript matrix with samples as observations. "
            "This is a qualitative sample-level view and should not be overinterpreted given the small sample size."
        )

        pca_df, explained = run_pca(expr)
        fig = px.scatter(
            pca_df,
            x="PC1",
            y="PC2",
            color="group",
            text="sample",
            title="Sample-level PCA (DMD vs WT)",
            color_discrete_map={"DMD": "#D1495B", "WT": "#2E86AB"},
        )
        fig.update_traces(textposition="top center", marker=dict(size=13, line=dict(width=1, color="white")))
        fig.update_layout(
            xaxis_title=f"PC1 ({explained[0]:.2f}% variance)",
            yaxis_title=f"PC2 ({explained[1]:.2f}% variance)",
        )
        st.plotly_chart(fig, use_container_width=True)

        st.info("Do not treat this PCA as inferential evidence; it is an exploratory visual summary.")

    elif section == "Figure Gallery":
        st.header("Figure Gallery")
        fig_dir = project_root / "downstream_analysis" / "figures"

        gallery = [
            (
                "Volcano-like Plot",
                fig_dir / "volcano_like_plot.png",
                "Descriptive transcript-level contrast using log2FC magnitude (no p-values/FDR).",
            ),
            (
                "Top Transcript Heatmap",
                fig_dir / "top_transcripts_heatmap.png",
                "Top |log2FC| transcripts with row-scaled log2(TPM + 1) across DMD/WT samples.",
            ),
            (
                "PCA Plot",
                fig_dir / "pca_plot.png",
                "Sample-level exploratory PCA from log2(TPM + 1) transformed matrix.",
            ),
            (
                "Summary Panel",
                fig_dir / "summary_panel.png",
                "Combined visual overview (volcano-like, PCA, and top transcript bars).",
            ),
        ]

        for title, path, caption in gallery:
            st.subheader(title)
            st.caption(caption)
            if path.exists():
                st.image(str(path), use_container_width=True)
            else:
                st.warning(f"Missing figure: {path}")

    elif section == "Reproducibility":
        st.header("Reproducibility")
        st.code("""pip install -r requirements.txt\nstreamlit run app.py""", language="bash")
        st.markdown(
            "Input files used by this app:\n"
            "- `results/matrix/expression_matrix.tsv`\n"
            "- `results/matrix/dmd_vs_wt_summary.tsv`"
        )

    elif section == "Limitations":
        st.header("Scientific Limitations")
        st.markdown(
            "- Small sample size (2 DMD, 2 WT).\n"
            "- TPM-based descriptive comparison only.\n"
            "- No p-values or FDR-adjusted significance estimates.\n"
            "- Transcript-level IDs are not gene-symbol annotated unless external mapping is provided.\n"
            "- Not for clinical interpretation."
        )


if __name__ == "__main__":
    main()
