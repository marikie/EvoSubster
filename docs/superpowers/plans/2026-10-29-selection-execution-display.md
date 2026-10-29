# Selection execution display

- [x] Capture the actual R invocation and every effective CLI value. Mark sources
  as specified, default, inherited from --idt-threshold, or derived cache path.
- [x] Replace the final tree's threshold heading with command, selection settings,
  runtime settings and paths, while retaining genus exemption and role labels.
- [x] Wrap long text using PDF font measurements and grow page height/top margin.
- [x] Test provenance, quoting, main-to-PDF wiring, long paths, missing metadata,
  PDF text/geometry, and existing selection/tree behavior. Inspect a rendered PDF.
- [x] Update README, review and validate the combined feature branch changes.

The current feature/newick-input tree implementation is the baseline. Selection
rules, TSVs, Stage 0 and Newick labels stay unchanged. Environment-only LAST
settings are outside scope. Existing PDFs are not regenerated automatically.
