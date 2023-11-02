context("Test of R function karlin")

test_that("Example from the R documentation of karlin(.)", {
  p <- karlin(150, 10000, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -5, 5)
  expect_equal(p, 0.03464975)
})

# Below : this check is running fine except on M1Mac configuration (polynomial root finding issue which throw a controlled error)
# test_that("Example from the a pig dataset (unit score scheme)", {
#   score1.prob.ext = c(`-2` = 0.0681561425858712, `-1` = 0.759474076388722, `0` = 0.121671514073504, 
#                       `1` = 0.0310891041511429, `2` = 0.0125039059683196, `3` = 0.00563901641708531,
#                       `4` = 0.0010912674566738, `5` = 0.000242770954017739, `6` = 7.93212027978752e-05,
#                       `7` = 3.60550921808523e-05, `8` = 1.44220368723409e-05, `9` = 2.40367281205682e-06)
#   p <- karlin(248, 416030, score1.prob.ext, -2, 9)
#   expect_equal(p, 0.0)
# })
