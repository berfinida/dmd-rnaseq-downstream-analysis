#!/usr/bin/env python3
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

SAMPLE_COLS = ["DMD1", "DMD2", "WT1", "WT2"]


def load_data(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr_path = project_root / "results" / "matrix" / "expression_matrix.tsv"
    summary_path = project_root / "results" / "matrix" / "dmd_vs_wt_summary.tsv"

    expr = pd.read_csv(expr_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")
    return expr, summary


def build_filtered_table(summary: pd.DataFrame, min_abs_log2fc: float, direction: str) -> pd.DataFrame:
    table_df = summary.copy()
    table_df["abs_log2FC"] = table_df["log2FC"].abs()
    table_df = table_df[table_df["abs_log2FC"] >= min_abs_log2fc]

    if direction == "upregulated":
        table_df = table_df[table_df["log2FC"] > 0].sort_values("log2FC", ascending=False)
    elif direction == "downregulated":
        table_df = table_df[table_df["log2FC"] < 0].sort_values("log2FC", ascending=True)
    else:
        table_df = table_df.sort_values("abs_log2FC", ascending=False)

    return table_df


def main() -> None:
    st.set_page_config(
        page_title="DMD vs WT Transcriptomics Explorer",
        page_icon="🧬",
        layout="wide",
    )

    project_root = Path(__file__).resolve().parent
    expr, summary = load_data(project_root)

    if "log2FC" not in summary.columns or "Name" not in summary.columns:
        st.error("Summary file must include 'Name' and 'log2FC' columns.")
        return

    st.title("DMD vs WT Transcriptomics Explorer")
    st.caption(
        "This app visualizes TPM-based descriptive summaries and does not perform formal statistical differential expression testing."
    )

    st.sidebar.title("Navigation")
    page = st.sidebar.radio(
        "Go to",
        [
            "Overview",
            "Transcript Explorer",
            "Top Transcript Tables",
            "Interactive Charts",
            "Figure Gallery",
        ],
    )

    st.sidebar.header("Filters")
    min_abs_log2fc = st.sidebar.slider(
        "Minimum absolute log2FC", min_value=0.0, max_value=10.0, value=2.0, step=0.1
    )
    top_n = st.sidebar.slider("Top N transcripts", min_value=5, max_value=100, value=20, step=5)
    direction = st.sidebar.selectbox(
        "Direction", options=["upregulated", "downregulated", "both"], index=2
    )

    filtered = build_filtered_table(summary, min_abs_log2fc, direction)

    if page == "Overview":
        st.header("Project Overview")
        st.write(
            "Interactive exploratory dashboard for TPM-based transcript summaries across "
            "WT1, WT2, DMD1, and DMD2 samples."
        )

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Number of transcripts", f"{len(summary):,}")
        col2.metric("Number of samples", str(sum(c in expr.columns for c in SAMPLE_COLS)))
        col3.metric("Max positive log2FC", f"{summary['log2FC'].max():.3f}")
        col4.metric("Max negative log2FC", f"{summary['log2FC'].min():.3f}")

        st.subheader("Expression Matrix Preview")
        st.dataframe(expr.head(25), use_container_width=True)

    elif page == "Transcript Explorer":
        st.header("Transcript Explorer")

        top_candidates = (
            summary.assign(abs_log2FC=summary["log2FC"].abs())
            .sort_values("abs_log2FC", ascending=False)
            .head(100)["Name"]
            .astype(str)
            .tolist()
        )

        transcript_query = st.text_input("Search by Ensembl transcript ID", placeholder="ENSMUST...")
        dropdown_choice = st.selectbox("Or choose from top transcripts", options=["(none)"] + top_candidates)

        selected_id = transcript_query.strip() if transcript_query.strip() else None
        if not selected_id and dropdown_choice != "(none)":
            selected_id = dropdown_choice

        if selected_id:
            expr_match = expr[expr["Name"].astype(str) == selected_id]
            summary_match = summary[summary["Name"].astype(str) == selected_id]

            if expr_match.empty:
                st.warning("Transcript ID not found.")
            else:
                st.subheader(f"Selected transcript: {selected_id}")
                tpm_vals = expr_match.iloc[0]
                available_samples = [c for c in SAMPLE_COLS if c in expr_match.columns]

                show_cols = ["Name", *available_samples]
                st.dataframe(expr_match[show_cols], use_container_width=True)

                if not summary_match.empty:
                    cols = [c for c in ["Name", "DMD_mean_TPM", "WT_mean_TPM", "log2FC"] if c in summary_match.columns]
                    st.dataframe(summary_match[cols], use_container_width=True)

                plot_df = pd.DataFrame(
                    {
                        "sample": available_samples,
                        "TPM": [float(tpm_vals[s]) for s in available_samples],
                    }
                )
                plot_df["group"] = plot_df["sample"].str.extract(r"(DMD|WT)")

                fig_line = px.line(
                    plot_df,
                    x="sample",
                    y="TPM",
                    markers=True,
                    color="group",
                    title="Selected transcript TPM across samples",
                    color_discrete_map={"DMD": "#D1495B", "WT": "#2E86AB"},
                )
                st.plotly_chart(fig_line, use_container_width=True)

                fig_bar = px.bar(
                    plot_df,
                    x="sample",
                    y="TPM",
                    color="group",
                    title="Selected transcript TPM bar chart",
                    color_discrete_map={"DMD": "#D1495B", "WT": "#2E86AB"},
                )
                st.plotly_chart(fig_bar, use_container_width=True)

    elif page == "Top Transcript Tables":
        st.header("Interactive Top Upregulated/Downregulated Table")
        st.dataframe(filtered.head(top_n), use_container_width=True)

    elif page == "Interactive Charts":
        st.header("Interactive Charts")

        volcano_df = summary.copy()
        volcano_df["abs_log2FC"] = volcano_df["log2FC"].abs()
        volcano_df["highlight"] = volcano_df["abs_log2FC"] >= min_abs_log2fc

        fig_volcano = px.scatter(
            volcano_df,
            x="log2FC",
            y="abs_log2FC",
            color="highlight",
            hover_data=["Name"],
            title="Volcano-like plot (interactive)",
            color_discrete_map={True: "#E45756", False: "#4C78A8"},
        )
        fig_volcano.update_traces(marker=dict(size=6, opacity=0.65))
        fig_volcano.add_vline(x=0, line_dash="dash", line_color="gray")
        st.plotly_chart(fig_volcano, use_container_width=True)

        up_df = summary[summary["log2FC"] > 0].sort_values("log2FC", ascending=False).head(top_n)
        down_df = summary[summary["log2FC"] < 0].sort_values("log2FC", ascending=True).head(top_n)

        c1, c2 = st.columns(2)
        with c1:
            fig_up = px.bar(
                up_df,
                x="log2FC",
                y="Name",
                orientation="h",
                title=f"Top {top_n} upregulated transcripts",
                color_discrete_sequence=["#D1495B"],
            )
            fig_up.update_layout(yaxis=dict(autorange="reversed"))
            st.plotly_chart(fig_up, use_container_width=True)

        with c2:
            fig_down = px.bar(
                down_df,
                x="log2FC",
                y="Name",
                orientation="h",
                title=f"Top {top_n} downregulated transcripts",
                color_discrete_sequence=["#2E86AB"],
            )
            fig_down.update_layout(yaxis=dict(autorange="reversed"))
            st.plotly_chart(fig_down, use_container_width=True)

    elif page == "Figure Gallery":
        st.header("Figure Gallery")
        fig_dir = project_root / "downstream_analysis" / "figures"

        figures = [
            ("Volcano-like Plot", fig_dir / "volcano_like_plot.png"),
            ("Top Transcript Heatmap", fig_dir / "top_transcripts_heatmap.png"),
            ("PCA Plot", fig_dir / "pca_plot.png"),
            ("Summary Panel", fig_dir / "summary_panel.png"),
        ]

        for title, fig_path in figures:
            st.subheader(title)
            if fig_path.exists():
                st.image(str(fig_path), use_container_width=True)
            else:
                st.info(f"Figure not found: {fig_path}")


if __name__ == "__main__":
    main()
