test_that("Reused subject IDs raises error.", {
  data <- data.frame(
    time = c(1, 1),
    status = c(1, 2),
    idx = c(1, 1),
    arm = c(1, 0)
  )
  
  expect_error(
    CompareAUCs(
      data$idx,
      data$time,
      data$status,
      data$arm,
      tau = 1
    )
  )
})

test_that("Test AUC comparison.", {
  
  tau <- 4
  
  GetAUCDiff <- function(data){
    out <-suppressWarnings(
      CompareAUCs(
        data$idx,
        data$time,
        data$status,
        data$arm,
        tau = tau,
        cens_after_last = FALSE
      )
    )
  }
  
  
  # Case 1.
  data1 <- data.frame(
    time = c(1, 1, 2, 2),
    status = c(1, 1, 1, 1),
    idx = c(1, 2, 3, 4)
  )
  data1$arm <- 1
  auc1 <- 0.5 * (2 - 1) + 1.0 * (tau - 2)
  
  data0 <- data.frame(
    time = c(1, 1, 2, 2),
    status = c(0, 2, 0, 2),
    idx = c(5, 6, 7, 8)
  )
  data0$arm <- 0
  auc0 <- 0
  
  data <- rbind(data1, data0)
  
  observed <- GetAUCDiff(data)
  expect_equal(observed@MargAreas$area[1], auc0)
  expect_equal(observed@MargAreas$area[2], auc1)
  expect_equal(observed@Pvals$observed[1], auc1 - auc0)
  
  
  # Case 2.
  data1 <- data.frame(
    time = c(1, 2, 3, 4),
    status =c(1, 2, 2, 1),
    idx = c(1, 2, 3, 4)
  )
  data1$arm <- 1
  auc1 <- 0.25 * (4 - 1) + 0.5 * (tau - 4)
  
  data0 <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(2, 2, 1, 0),
    idx = c(5, 6, 7, 8)
  )
  data0$arm <- 0
  auc0 <- 0.25 * (tau - 3)
  
  data <- rbind(data1, data0)
  
  observed <- GetAUCDiff(data)
  expect_equal(observed@MargAreas$area[1], auc0)
  expect_equal(observed@MargAreas$area[2], auc1)
  expect_equal(observed@Pvals$observed[1], auc1 - auc0)

})