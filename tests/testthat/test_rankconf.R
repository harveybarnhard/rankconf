library(rankconf)

df = read.csv(
  "https://opportunityinsights.org/wp-content/uploads/2018/10/cz_outcomes_simple.csv"
)
cols = c("cz", "czname", "kfr_pooled_pooled_p25", "kfr_pooled_pooled_p25_se")
df = df[1:20, cols]
n  = nrow(df)
colnames(df) = c("cz", "czname", "kfr_p25", "kfr_p25_se")

typel = list(
  c("PCER", "NAIVE"),
  c("FWER", "B"),
  c("FWER", "HB"),
  c("FWER", "R"),
  c("FDR", "BY")
)

for(i in 1:length(typel)) {
  test_that(paste0("Testing type: ", typel[[i]][1], " with method: ", typel[[i]][2]), {
    output = rankconf(
      df$kfr_p25, sig2=df$kfr_p25_se^2, type=typel[[i]][1], method=typel[[i]][2]
    )
    expect_true(all(output$L <= output$U))
    expect_true(all(output$L >= 1 & output$L <= n))
    expect_true(all(output$U >= 1 & output$U <= n))
    expect_false(anyNA(output$L))
    expect_false(anyNA(output$U))
    expect_length(output$L, n)
    expect_length(output$U, n)
  })
}
