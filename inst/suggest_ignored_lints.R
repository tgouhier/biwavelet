# Run this script if you want to ignore lints that currently appear in your code

library(magrittr)
library(dplyr)

lintr::lint_package() %>%
as.data.frame %>%
group_by(linter) %>%
tally(sort = TRUE) %$%
sprintf("linters: with_defaults(\n    %s\n    NULL\n  )\n",
        paste0(linter, " = NULL, # ", n, collapse = "\n    ")) %>%
cat()
