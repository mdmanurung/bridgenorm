# Purpose: Build the Quarto documentation site with altdoc into docs/, then fail if private material leaked into
#          the site or the public sources.
# Run:     Rscript tools/build-site.R   (from the package root; needs quarto on PATH and bridgenorm installed)
# Inputs:
#   DESCRIPTION, README.md, NEWS.md, man/, vignettes/    package sources
#   altdoc/quarto_website.yml                             navbar, theme and sidebar settings
# Outputs:
#   docs/        the site (gitignored on main, deployed from gh-pages); open docs/index.html
#   _quarto/     altdoc's staging copy, removed at the end

# Patterns that must never appear in public files: private notes and folders, local paths, real marker names
# (CD followed by digits), and decision ids of the private log. Extra private terms (data-set names, ID
# formats) live in an untracked file so this script does not publish them.
PRIVATE <- c(
  "HANDOFF",
  "experiments/",
  "plans/",
  "/exports/",
  "\\bCD[0-9]+",
  "\\b(LG|ZS|NS|RV|VC|RL)[0-9]+[a-z]?\\b"
)
PRIVATE_FILE <- "experiments/leak-patterns.txt" # one extended regex per line; absent on other machines
PUBLIC <- c(
  "docs",
  "README.md",
  "NEWS.md",
  "man",
  "vignettes",
  "tests",
  "DESCRIPTION"
)

# Build ------------------------------------------------------------------------

# freeze = FALSE reruns every chunk, so figures always match the code.
altdoc::render_docs(freeze = FALSE, verbose = TRUE)
unlink("_quarto", recursive = TRUE)

# Check for private material ------------------------------------------------

patterns <- c(PRIVATE, if (file.exists(PRIVATE_FILE)) readLines(PRIVATE_FILE))
patterns <- patterns[nzchar(patterns)]
leaks <- unlist(lapply(patterns, function(p) {
  hits <- suppressWarnings(system2(
    "grep",
    c("-rliE", shQuote(p), PUBLIC),
    stdout = TRUE
  )) # exit 1 = no match
  if (length(hits) > 0) paste0(p, ": ", paste(hits, collapse = ", "))
}))
if (length(leaks) > 0) {
  stop("Private material in public files:\n", paste(leaks, collapse = "\n"))
}
