# Selection identity thresholds

Approved specification: the conversation's revised plan, 2026-10-28.

- [x] Extend the shared Python filter with independent ingroup/outgroup thresholds
  (80 each) and `--min-idt-gap-between-out-in` (0). Both gaps must strictly
  exceed the configured value unless the existing target genus pattern exempts
  them. Absolute thresholds and the two-vs-one exclusion remain mandatory.
- [x] Wire the options through the R selector and tree-mode shell wrapper.
  Preserve `--idt-threshold` as a common fallback; explicit individual values
  take precedence regardless of argument order. Validate finite values in [0,100].
- [x] Add applied thresholds, both gaps, gap verdict and exemption to verdict TSVs;
  document CLI behavior and verify Python, R, wrapper and report regressions.

Report CLI expansion is outside scope. Existing Python report calls retain the
common threshold and zero-gap defaults. Work in the existing isolated worktree.
