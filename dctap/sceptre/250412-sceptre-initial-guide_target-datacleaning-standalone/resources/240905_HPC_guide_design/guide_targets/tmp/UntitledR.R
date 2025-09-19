library(dplyr)
df <- data.frame(
  grna_target = discovery_pairs$grna_target,
  response_id = discovery_pairs$response_id,
  stringsAsFactors = FALSE
)

write.table(df, file = here("tmp/discovery_pairs.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

df1 <- data.frame(
  grna_target = positive_pairs$grna_target,
  response_id = positive_pairs$response_id,
  stringsAsFactors = FALSE
)

write.table(df1, file = here("tmp/positive_pairs.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)