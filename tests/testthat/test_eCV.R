
test_that("eCV function works as expected", {
    # Test with valid input
    x <- matrix(c(1, 2, 3, 4), nrow = 2)
    result <- eCV(x)
    expect_equal(length(result), 2)
    expect_true(all(names(result) %in% c("eCV", "post_prob")))

    # Test with invalid input
    expect_error(eCV(x = "invalid input"))
    expect_error(eCV(x, max.ite = "invalid"))
    expect_error(eCV(x, n_threads = 0))
})
