library(testthat)
library(supervisoR)

test_that("sanitize_filename works correctly", {
  expect_equal(sanitize_filename("path/to/file name with spaces.txt"), "path_to_file_name_with_spaces_txt")
  expect_equal(sanitize_filename("filename"), "filename")
  expect_equal(sanitize_filename("file@name#1!.txt"), "file_name_1__txt")
})

test_that("is_path correctly identifies paths", {
  expect_true(is_path("path/to/file.txt"))
  expect_true(is_path("C:\\Users\\User\\file.txt"))
  expect_false(is_path("filename"))
  expect_false(is_path("file.name"))
})
