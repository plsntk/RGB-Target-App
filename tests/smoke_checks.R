Sys.setenv(RGB_TARGET_APP_DISABLE_AUTORUN = "1")
source("app.R", local = TRUE)

# compute_grid_paged: should detect overflow when pages exceed max_pages
res_overflow <- compute_grid_paged(
  page_w_mm = 297,
  page_h_mm = 420,
  patch_mm = 20,
  gap_mm = 2,
  margin_side_mm = 4,
  margin_leading_mm = 23,
  margin_trailing_mm = 33,
  max_chart_w_mm = 300,
  n_patches = 5000,
  max_pages = 1
)
stopifnot(isFALSE(res_overflow$fits))

# pad_to_capacity: exact capacity + alternating fill
base_dt <- data.table::data.table(R = c(1L, 2L), G = c(1L, 2L), B = c(1L, 2L))
padded <- pad_to_capacity(base_dt, capacity = 6, mode = "alternate")
stopifnot(nrow(padded) == 6)
stopifnot(identical(as.integer(padded$R[3:6]), c(255L, 0L, 255L, 0L)))

# make_offgrey_rings: honors override and preserves bounds
neutrals <- make_neutrals(5)
rings <- make_offgrey_rings(neutrals, rings = 2, total_colors = 6000, delta_max_override = 4)
stopifnot(identical(attr(rings, "delta_max_used"), 4L))
stopifnot(all(rings$R >= 0L & rings$R <= 255L))
stopifnot(all(rings$G >= 0L & rings$G <= 255L))
stopifnot(all(rings$B >= 0L & rings$B <= 255L))

# deterministic fill selection for same seed
f1 <- make_fill_poisson_oklab_balanced(50, seed = 123)
f2 <- make_fill_poisson_oklab_balanced(50, seed = 123)
stopifnot(identical(f1, f2))

cat("All smoke checks passed.\n")
