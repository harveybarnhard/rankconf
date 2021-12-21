library(rankconf)

df = read.csv(
  "https://opportunityinsights.org/wp-content/uploads/2018/10/cz_outcomes_simple.csv"
)
cols = c("cz", "czname", "kfr_pooled_pooled_p25", "kfr_pooled_pooled_p25_se")
df = df[, cols]
n  = nrow(df)

test_that("Testing naive method", {
  output = rankconf(df$kfr_p25, sig2=df$kfr_p25_se^2, type="PCER", method="NAIVE")
  expect_true(all(output$L <= output$U))
  expect_false(anyNA(output$L))
  expect_false(anyNA(output$U))
  expect_length(length(output$L), n)
  expect_length(length(output$U), n)
})
