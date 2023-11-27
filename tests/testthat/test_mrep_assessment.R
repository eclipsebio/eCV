test_that("Testing if mIDR and gIDR are the same as IDR
 for n_reps = 2", {
    set.seed(42)
    library("eCV")
    library("testthat")
    out1 <- simulate_data(scenario = 1, n_features = 500)

    params <- list(
        eCV = list(max.ite = 100),
        IDR = list(
            mu = 2,
            sigma = 1.3,
            rho = 0.8,
            p = 0.7,
            eps = 1e-3,
            max.ite = 50
        ),
        gIDR = list(
            mu = 2,
            sigma = 1.3,
            rho = 0.8,
            p = 0.7,
            eps = 1e-3,
            max.ite = 50
        ),
        mIDR = list(
            mu = 2,
            sigma = 1.3,
            rho = 0.8,
            p = 0.7,
            eps = 1e-3,
            max.ite = 50
        )
    )

    idr_res <- mrep_assessment(
        x = out1$sim_data,
        method = "IDR",
        param = params[["IDR"]]
    )$rep_index

    gidr_res <- mrep_assessment(
        x = out1$sim_data,
        method = "gIDR",
        param = params[["gIDR"]]
    )$rep_index

    midr_res <- mrep_assessment(
        x = out1$sim_data,
        method = "mIDR",
        param = params[["mIDR"]]
    )$rep_index
    expect_equal(idr_res, gidr_res)
    expect_equal(idr_res, midr_res)
})
