# Vignette
usethis::use_vignette("EMOGEA")
devtools::build_vignettes()

# Build doc and NAMESPACE
golem::detach_all_attached()
devtools::document()

# Check code coverage and see what still needs to be tested
golem::detach_all_attached()
devtools::test()
golem::detach_all_attached()
devtools::test_coverage()

# Check the package
golem::detach_all_attached()
devtools::check()

# Test if you can install the package
golem::detach_all_attached()
remotes::install_local(upgrade = "never", force = TRUE, build_vignettes = TRUE)

# Create a tar.gz
devtools::build()
