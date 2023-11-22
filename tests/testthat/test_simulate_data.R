
test_that("Testing if
          simulations work properly ", {
    set.seed(42)
    library("eCV")
    out1 <- simulate_data(scenario = 1)
    expect_equal(out1$sim_params$feature_group, rep(c(1, 2), c(3500, 6500)))
    out2 <- simulate_data(scenario = 2)
    expect_equal(out2$sim_params$feature_group, rep(c(1, 2), c(7000, 3000 )))
    out3 <- simulate_data(scenario = 3)
    expect_equal(out3$sim_params$feature_group, rep(c(1, 2), c(9500, 500 )))
    out4 <- simulate_data(scenario = 4)
    expect_equal(out4$sim_params$feature_group, rep(c(1, 2, 3), c(2800, 6500, 700 )))
})
