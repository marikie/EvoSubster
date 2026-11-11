"""Verify rendered alpha and preserve white data fills across figure producers."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("Rscript") and shutil.which("pdftocairo"),
                     "Rendered transparency checks require Rscript and pdftocairo")
class TransparentFigureTests(unittest.TestCase):
    def test_rendered_producers(self):
        work = Path(tempfile.mkdtemp(prefix="transparent-render-", dir=ROOT / "tmp/codex"))
        env = dict(os.environ, EVOSUBSTER_FONT_CACHE=str(ROOT / "tmp/codex-checks/font-cache"))

        def run(command):
            result = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + "\n" + result.stderr)

        def check_alpha(path):
            with Image.open(path) as image:
                # Low-colour PNGs may store alpha in their palette's tRNS chunk.
                alpha = image.convert("RGBA").getchannel("A")
                self.assertEqual(alpha.getextrema()[0], 0, str(path))
                if path.name not in ("empty_metrics.png", "Rplots_raster.png"):
                    self.assertGreater(alpha.getextrema()[1], 0, str(path))
                self.assertEqual(alpha.getpixel((0, 0)), 0, str(path))
                if path.name in ("metrics.png", "scatter.png"):
                    self.assertGreater(alpha.histogram()[0] / (image.width * image.height), 0.5,
                                       "The plot panels must be transparent too")

        code = r'''
args <- commandArgs(TRUE); root <- args[1]; out <- args[2]
source(file.path(root, "src/visualize/pca_2d.R"))
original_theme <- ggplot2::theme_get()
white <- ggplot2::ggplot(data.frame(x = 1:2, y = 1:2), ggplot2::aes(x, y)) +
  ggplot2::geom_col(fill = "white", colour = "black") + ggplot2::theme_bw() + transparent_background_theme()
ggplot2::ggsave(file.path(out, "white_data.png"), white, width = 4, height = 3, bg = "transparent")
ggplot2::ggsave(file.path(out, "white_data.pdf"), white, width = 4, height = 3, bg = "transparent")
stopifnot(identical(ggplot2::theme_get(), original_theme))
plot_metrics(file.path(out, "metrics.png"), c(10, 5, 3), c(0.2, 0.5, 0.3))
plot_metrics(file.path(out, "empty_metrics.png"), numeric(), numeric())
write_placeholder_plot(file.path(out, "placeholder.png"), "No samples")
coordinates <- matrix(c(0,0,2,0,0,2,2,2,10,10,12,10,10,12,12,12), ncol = 2, byrow = TRUE)
plot_pca_scatter(file.path(out, "scatter.png"), coordinates, letters[1:8],
  rep(c("Mucoromycota", "Apicomplexa"), each = 4), rep(1:2, each = 4),
  "phylum", "test", cluster_boundary = "hull", show_labels = FALSE)
source(file.path(root, "src/visualize/spectrum_shape.R"))
x <- sin(outer(1:8, 1:96, "*") / 10); x[, 1] <- 0
rownames(x) <- paste0("S", 1:8); colnames(x) <- paste0("channel", 1:96)
meta <- data.frame(sample_id = rownames(x), species = paste("Species", 1:8), display_id = rownames(x),
  phylum = rep(c("Mucoromycota", "Apicomplexa"), each = 4))
shape_plots(x, meta, shape_cluster(x, "ward", "2", 4), coordinates, "phylum", out)
source(file.path(root, "src/select/trio_tree_outputs.R"))
tree <- ape::read.tree(text = "(A:1,(B:1,C:1):1);")
write_stage0_tree(tree, data.frame(tip = c("A","B","C"), accession = c("A1","B2","C3")), out)
lines <- readLines(file.path(root, "src/report/report_template.Rmd"))
eval(parse(text = lines[grep("^render_pdf_preview <-", lines):(grep("^value_or_na <-", lines) - 1)]))
has_pdftools <- TRUE
cache <- file.path(out, "preview_cache"); images <- file.path(out, "preview_images")
dir.create(cache); dir.create(images)
pdf <- file.path(out, "white_data.pdf")
writeLines("stale flattened cache", file.path(cache, paste0(unname(tools::md5sum(pdf)), "-1.png")))
preview <- render_pdf_preview(pdf, cache, images, dpi = 36)
stopifnot(file.exists(preview), grepl("-transparent-v1.png$", preview))
bench <- data.frame(trio_label = paste0("T", 1:6), lineage = rep(c("fungi", "cnidaria"), each = 3),
  wall_clock = paste0(3:8, ":00"), cpu_percent = 800, max_rss_kb = 1000000)
metadata <- data.frame(trio_label = bench$trio_label, lineage = bench$lineage,
  size_org1_bp = (2:7)*10000000, size_org2_bp = 20000000, size_org3_bp = 20000000,
  idt_12 = 85:90, idt_13 = 86:91, idt_23 = 97)
write.table(bench, file.path(out, "bench.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(metadata, file.path(out, "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
'''
        run(["Rscript", "-e", code, str(ROOT), str(work)])
        run(["Rscript", str(ROOT / "src/visualize/bench_figure.R"), str(work / "bench.tsv"), str(work / "metadata.tsv"), str(work)])
        fixture = ROOT / "test/fixtures/pipeline_counts/20261031101732/Oikalb2"
        single = next(p for p in (fixture / "singlenuc").glob("*.tsv") if "_ncds" not in p.name)
        dinuc = next(p for p in (fixture / "dinuc").glob("*.tsv") if "_ncds" not in p.name)
        for script, source in (("sbmut.R", single), ("logRatioPlot.R", single), ("dinucleotide-plot.R", dinuc)):
            run(["Rscript", str(ROOT / "src/visualize" / script), str(source), str(work), "Test species"])
        for path in work.rglob("*.png"):
            if not path.name.endswith("-1.png"):
                check_alpha(path)
        for path in work.rglob("*.pdf"):
            prefix = work / (path.stem + "_raster")
            run(["pdftocairo", "-png", "-transp", "-singlefile", "-f", "1", "-l", "1", "-r", "36", str(path), str(prefix)])
            check_alpha(Path(str(prefix) + ".png"))
        for name, minimum in (("white_data.png", 1000), ("heatmap.png", 100)):
            with Image.open(work / name) as image:
                colors = image.convert("RGBA").getcolors(image.width * image.height)
                count = sum(n for n, color in colors if color == (255, 255, 255, 255))
                self.assertGreater(count, minimum, f"White data fills must remain opaque: {name}")
        shutil.rmtree(work)  # Failed tasks retain their render evidence above.
