# Purpose: Build the local pkgdown site into docs/, then drop pages pkgdown makes from private root notes.
# Run:     Rscript pkgdown/build-site.R   (from the package root, R4_51)
# Inputs:
#   DESCRIPTION, README.md, NEWS.md, man/, vignettes/, _pkgdown.yml
# Outputs:
#   docs/   the site (gitignored); open docs/index.html

PRIVATE <- c("HANDOFF") # root .md files pkgdown renders although .Rbuildignore lists them

# Build ------------------------------------------------------------------------

pkgdown::clean_site(force = TRUE) # drop pages of functions that no longer exist
pkgdown::build_site(preview = FALSE, new_process = TRUE)

# Drop private pages -----------------------------------------------------------

# Remove the rendered pages, their search entries and their sitemap rows.
pattern <- paste(PRIVATE, collapse = "|")
file.remove(list.files("docs", pattern = paste0("^(", pattern, ")\\.(html|md)$"), full.names = TRUE))
search <- jsonlite::read_json("docs/search.json")
kept <- Filter(\(e) !grepl(pattern, jsonlite::toJSON(e, auto_unbox = TRUE)), search)
jsonlite::write_json(kept, "docs/search.json", auto_unbox = TRUE)
sitemap <- readLines("docs/sitemap.xml")
writeLines(sitemap[!grepl(pattern, sitemap)], "docs/sitemap.xml")
stopifnot("no private page left in docs" = length(system2("grep", c("-rl", pattern, "docs"), stdout = TRUE)) == 0)
