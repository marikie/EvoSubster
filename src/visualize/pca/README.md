`../pca_2d.R` remains the CLI and sourceable entry point. It loads these modules
into the caller's environment; invoke the entry point rather than sourcing
individual modules.

- `config.R`: dependencies, constants, and CLI parsing.
- `input.R`: run discovery, manifest/TSV loading, and ordered sample collection.
- `taxonomy.R`: metadata and taxonomy labels.
- `identity.R`: the shared Python identity-filter bridge.
- `compute.R`: observed/logratio features, scaling, PCA, and cluster selection.
- `plot.R`: standard scatter plots, boundaries, and metric plots.
- `output.R`: result tables, loadings, and run/filter audit persistence.
- `pipeline.R`: execution order and coordination of the modules.

Sample collection updates a shared audit environment as it proceeds. The entry
point's `main()` registers an exit handler before collecting inputs, preserving
completed trio decisions even when a later input fails. Computation and plot
calls retain their original order, including random seeds.
