# PPI Analysis

This repository provides a comprehensive pipeline for analyzing protein-protein interaction (PPI) networks, gene expression, domain enrichment, and biomarker identification across biological subtypes. The project is organized for modular, reproducible analysis in R, with clear separation of data, functions, scripts, and visualization outputs.

---

## Repository Structure

```
analysis/
  ├── full_data_ppi.R
  ├── data_gathering/
  │     ├── 01_get_GO_genes.R
  │     ├── 02_data_subsets.R
  │     └── 03_subtype_biomarkers.R
  ├── enrichment/
  │     └── enrichment_analysis.R
  ├── ppi_hub/
  │     ├── 04_ppi_network.R
  │     ├── 05_subtype_hubs.R
  │     ├── direct_overlap.R
  │     ├── hub_viz.R
  │     └── ppi_network_viz.R
  └── tmp/
        └── domains/
              ├── 01_Extract_Domains_from_Hub_Genes.R
              ├── 02_Domain_Enrichment_Analysis.R
              ├── 03_Visualize_Domain_Architecture.R
              ├── 04_Subtype_Specific_Domain_Analysis.R
              ├── 5_Domain_Network_Analysis.R
              ├── cell_cycle_ppi.R
              ├── comp_hub.R
              └── plot_domains.R

data/
  ├── subtypes_samples.Rdata
  ├── bio_markers/
  ├── domains/
  ├── expr/
  ├── GO/
  ├── hubs/
  ├── ppi_network/
  └── predicted/

functions/
  ├── get_ppi.R
  ├── get_subtype_biomarkers.R
  ├── get_subtype_hubs.R
  ├── helper_get_expr_metrics.R
  └── plot_ppi_network.R

viz/
  ├── cell_cycle/
  └── luminal_diff/
```

---

## Main Components

### 1. `analysis/`
Contains all R scripts for data processing, enrichment analysis, PPI network construction, hub identification, and domain analysis.

- **data_gathering/**: Scripts for collecting gene ontology (GO) genes, creating data subsets, and identifying subtype biomarkers.
- **enrichment/**: Performs enrichment analysis on gene sets.
- **ppi_hub/**: Builds PPI networks, identifies hub genes, visualizes networks, and analyzes hub overlaps.
- **tmp/domains/**: Extracts domains from hub genes, performs domain enrichment, visualizes domain architectures, and analyzes subtype-specific domains.

### 2. `data/`
Stores all input and intermediate data files in RData format, organized by analysis type:
- **bio_markers/**: Biomarker data for subtypes.
- **domains/**: Domain data for subtypes.
- **expr/**: Expression matrices for various conditions.
- **GO/**: GO term data for cell cycle and luminal differentiation.
- **hubs/**: Hub gene matrices and comparisons.
- **ppi_network/**: PPI network data for subtypes.
- **predicted/**: Predicted subtype and sample data.

### 3. `functions/`
Reusable R functions for core analysis tasks:
- `get_ppi.R`: Functions for PPI network construction.
- `get_subtype_biomarkers.R`: Functions for biomarker identification.
- `get_subtype_hubs.R`: Functions for hub gene analysis.
- `helper_get_expr_metrics.R`: Helper functions for expression metrics.
- `plot_ppi_network.R`: Functions for network visualization.

### 4. `viz/`
Contains all output visualizations (PDFs) for cell cycle and luminal differentiation analyses, including:
- Combined PPI networks
- Domain architecture plots
- GO enrichment results
- Hub centrality vs. expression plots
- Hub specificity plots

---

## Usage

1. **Data Preparation**: Place all required `.Rdata` files in the appropriate subfolders under `data/`.
2. **Run Analysis Scripts**: Execute scripts in `analysis/` as needed, following the logical order (data gathering → enrichment → PPI network → hub/domain analysis).
3. **Use Functions**: Import functions from `functions/` in your scripts for modular analysis.
4. **View Results**: Output visualizations are saved in `viz/` for interpretation and publication.

---

## Typical Workflow

1. **Gene and Data Collection**: Use scripts in `analysis/data_gathering/` to collect GO genes and create data subsets.
2. **Biomarker Identification**: Identify subtype-specific biomarkers.
3. **Enrichment Analysis**: Run enrichment analysis on gene sets.
4. **PPI Network Construction**: Build PPI networks and identify hub genes.
5. **Domain Analysis**: Extract and analyze protein domains from hub genes.
6. **Visualization**: Generate publication-ready figures in `viz/`.

---

## Requirements

- R (version 4.x recommended)
- Required R packages (see individual scripts for dependencies; typically includes `igraph`, `ggplot2`, `dplyr`, etc.)

---

## Customization

- Add new data to `data/` and update scripts in `analysis/` as needed.
- Extend or modify functions in `functions/` for new analysis types.
- Save new visualizations in `viz/` for additional subtypes or conditions.

---

## Contact

For questions or contributions, please contact the repository owner.

---

## License

Specify your license here (e.g., MIT, GPL-3.0, etc.).

---

## Citation

If you use this repository in your research, please cite appropriately.

---

This README provides a high-level overview. For details on each analysis step, refer to the comments and documentation within individual R scripts.
