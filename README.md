Usage of the PseudoMeta R Package
The PseudoMeta package is designed to analyze single-cell RNA sequencing (scRNA-seq) data by estimating pseudotime-associated gene expression trends using non-parametric statistical methods.

Step 1: Create the PseudoMeta Object
Function: CreatePseudometa()

Input: A Seurat object (for scRNA-seq expression and dimensionality reduction) and a CDS object (from monocle3, used to extract pseudotime).

Process:

Selects the top Varnum variable genes.

Extracts the expression matrix for these genes.

Stores UMAP and (optionally) PCA reductions.

Retrieves pseudotime from the CDS object.

Output: A new PseudoMeta S4 object, with slots for:

assays: expression data and selected features

reductions: UMAP/PCA

meta.data: metadata including pseudotime

analysis: initialized list for future test results

Step 2: Calculate Gene Expression Trends Along Pseudotime
Function: Pseudotrend()

Input: A PseudoMeta object, statistical method name (Spearman, Kendall, Mann_Kendall, or Page), and a minimum number of expressing cells.

Process:

Filters cells with non-zero pseudotime.

Sorts cells by pseudotime.

Applies the chosen non-parametric trend test to each gene.

Statistical Methods Implemented
1. Spearman Rank Correlation
Function: spearman_test()

Calculates Spearman's rho between pseudotime order and gene expression ranks.

Indicates monotonicity direction: "positive" or "negative".

2. Kendall's Tau Correlation
Function: kendall_test()

Computes Kendall’s tau, measuring the strength and direction of association between time and gene expression ranks.

Results include tau-based p.value and monotonicity.

3. Mann-Kendall Trend Test
Function: mann_kendall_test()

A non-parametric test for monotonic trends.

Outputs:

tau (Kendall’s tau)

S (trend score)

varS (variance of S)

p.value, direction of trend.

4. Page Test (Page's L Test)
Function: page_test()

Detects ordered alternatives in ranked data (increasing or decreasing trends).

Computes:

L: Page’s test statistic

Z: normalized test score

p.value, direction of trend

Class Definition
PseudoMeta S4 Class:

r
slots = c(
  assays = "list",
  analysis = "list",
  reductions = "list",
  meta.data = "data.frame",
  ID = "character"
)
Summary
The PseudoMeta package facilitates pseudotime trend analysis in single-cell data by:

Integrating Seurat and monocle3 outputs.

Extracting highly variable genes.

Applying robust, non-parametric statistical tests to detect monotonic gene expression changes across pseudotime.

You can choose among Spearman, Kendall, Mann-Kendall, or Page tests depending on your needs for correlation or trend detection robustness.
