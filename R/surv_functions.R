library(tidyverse)
library(pillar)
# Input data:
# table:
my_data <- tibble(
  ID = c(1, 2, 3),
  Name = c("John", "Alice", "Bob"),
  Age = c(25, 30, 22)
)

# Add descriptions to columns
pillar::pillar_options(tibble::pillar_options("pillar.comment" = list(
  ID = "Unique identifier",
  Name = "Person's name",
  Age = "Person's age"
)))

# Print the tibble to see the descriptions
print(my_data)