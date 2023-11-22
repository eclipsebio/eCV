test_that("gIDR function works as expected", {
    set.seed(42)
    library("eCV")
    library("testthat")
    out <- simulate_data(scenario = 1, n_features = 500, n_reps = 3)
    mu <- 2
    sigma <- 1
    rho <- 0.5
    p <- 0.8
    result <- gIDR(out$sim_data, mu, sigma, rho, p)
    expect_equal(length(result), 3) # Check if the function returns a list with three elements
    expect_true(all(names(result) %in% c("est_param", "idr", "IDR"))) # Check if the list contains the expected elements

    # Add more test cases here to cover different scenarios

    # Test with invalid input
    expect_error(gIDR(x = out, mu, sigma, rho, p))
    expect_error(gIDR(x = out$sim_data, "invalid mu", sigma, rho, p))
    expect_error(gIDR(x = out$sim_data, mu, "invalid sigma", rho, p))
    expect_error(gIDR(x = out$sim_data, mu, sigma, "invalid rho", p))
    expect_error(gIDR(x = out$sim_data, mu, sigma, rho, "invalid p"))
})
