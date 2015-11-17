if (requireNamespace("lintr", quietly = TRUE)) {
  context("lintr static code analysis")
  test_that("Package Style", {

    dl <- lintr::default_linters

    # the following linters will be removed
    dl["spaces_left_parentheses_linter"] <- NULL

    lintr::expect_lint_free(linters = dl)
  })
}
