# Purpose: Build the Quarto documentation site with altdoc into docs/, then fail if a private note leaked in.
# Run:     Rscript tools/build-site.R   (from the package root; needs quarto on PATH and bridgenorm installed)
# Inputs:
#   DESCRIPTION, README.md, NEWS.md, man/, vignettes/    package sources
#   altdoc/quarto_website.yml                             navbar, theme and sidebar settings
# Outputs:
#   docs/        the site (gitignored on main, deployed from gh-pages); open docs/index.html
#   _quarto/     altdoc's staging copy, removed at the end

PRIVATE <- "HANDOFF" # root note that must never appear in the site

# Build ------------------------------------------------------------------------

# freeze = FALSE reruns every chunk, so figures always match the code.
altdoc::render_docs(freeze = FALSE, verbose = TRUE)
unlink("_quarto", recursive = TRUE)

# Check for private notes ------------------------------------------------------

leaks <- suppressWarnings(system2("grep", c("-rli", PRIVATE, "docs"), stdout = TRUE)) # grep exits 1 when nothing matches
if (length(leaks) > 0) stop("docs/ mentions ", PRIVATE, " in: ", paste(leaks, collapse = ", "))
