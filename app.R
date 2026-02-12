# High-quality RGB printer target generator for MYIROtools (XML + multi-page TIFF + CGATS)
# - Mandatory-first (neutrals, primaries/secondaries, ramps, off-grey rings)
# - Gamma 2.2 perceptual-ish spacing for steps
# - Fill via greedy Poisson-disk / maximin sampling in OKLab (approx, deterministic)
# - Keep overshoot if it still fits within Max Pages
# - Scramble print order (seeded), add repeated white/black, pad to full rows-per-page

suppressPackageStartupMessages({
  library(shiny)
  library(data.table)
  library(tiff)
  library(zip)
})

# -------------------- utilities --------------------

clip255 <- function(x) {
  d <- dim(x)
  y <- as.integer(round(x))
  y[y < 0L] <- 0L
  y[y > 255L] <- 255L
  if (!is.null(d)) dim(y) <- d
  y
}

unique_rgb_dt <- function(dt) {
  dt <- as.data.table(dt)
  setnames(dt, c("R","G","B"))
  dt[, `:=`(R=clip255(R), G=clip255(G), B=clip255(B))]
  unique(dt)
}

add_siblings <- function(dt) {
  a <- copy(dt)
  b <- dt[, .(R=G, G=B, B=R)]
  c <- dt[, .(R=B, G=R, B=G)]
  unique(rbindlist(list(a,b,c)))
}

# Gamma-spaced steps including endpoints, unique, length may be <= k if rounding collides (acceptable)
gamma_steps <- function(k, gamma=2.2) {
  if (k <= 1) return(0L)
  t <- seq(0, 1, length.out = k)
  v <- round(255 * (t^(1/gamma)))
  clip255(unique(v))
}

# Linear-light steps mapped back to sRGB code values (prevents highlight bias)
linear_steps_srgb <- function(k) {
  if (k <= 1) return(0L)
  t <- seq(0, 1, length.out = k)          # linear light
  srgb01 <- ifelse(t <= 0.0031308, 12.92*t, 1.055*(t^(1/2.4)) - 0.055)
  clip255(unique(as.integer(round(255 * srgb01))))
}

# -------------------- OKLab conversion --------------------
# sRGB (0..255) -> linear -> OKLab (per Björn Ottosson)
srgb_to_linear <- function(u01) {
  # u01 numeric in [0,1]
  ifelse(u01 <= 0.04045, u01/12.92, ((u01 + 0.055)/1.055)^2.4)
}

rgb255_to_oklab <- function(rgb) {
  # rgb: Nx3 integer 0..255 matrix (or length-3 vector)
  rgb <- as.matrix(rgb)
  if (is.null(dim(rgb))) rgb <- matrix(rgb, ncol = 3, byrow = TRUE)
  if (ncol(rgb) != 3) stop("rgb255_to_oklab expects Nx3 input")

  rgb01 <- rgb / 255
  # clamp without dropping dimensions
  rgb01[rgb01 < 0] <- 0
  rgb01[rgb01 > 1] <- 1

  r <- srgb_to_linear(rgb01[,1])
  g <- srgb_to_linear(rgb01[,2])
  b <- srgb_to_linear(rgb01[,3])

  # linear sRGB -> LMS
  l <- 0.4122214708*r + 0.5363325363*g + 0.0514459929*b
  m <- 0.2119034982*r + 0.6806995451*g + 0.1073969566*b
  s <- 0.0883024619*r + 0.2817188376*g + 0.6299787005*b

  l_ <- l^(1/3)
  m_ <- m^(1/3)
  s_ <- s^(1/3)

  L <- 0.2104542553*l_ + 0.7936177850*m_ - 0.0040720468*s_
  A <- 1.9779984951*l_ - 2.4285922050*m_ + 0.4505937099*s_
  B <- 0.0259040371*l_ + 0.7827717662*m_ - 0.8086757660*s_

  cbind(L, A, B)
}

# -------------------- mandatory set --------------------

make_neutrals <- function(k) {
  v <- linear_steps_srgb(k)
  data.table(R=v, G=v, B=v)
}

make_hue_ramps <- function(k) {
  # ramps with k-ish steps (unique rounding allowed)
  v <- linear_steps_srgb(k) / 255 / 255  # 0..1
  hues <- list(
    R = c(255,   0,   0),
    G = c(  0, 255,   0),
    B = c(  0,   0, 255),
    C = c(  0, 255, 255),
    M = c(255,   0, 255),
    Y = c(255, 255,   0)
  )

  out <- list()
  for (h in hues) {
    h <- as.numeric(h)

    # shade: black -> hue
    shade <- cbind(h[1]*v, h[2]*v, h[3]*v)

    # tint: white -> hue
    tint <- cbind(
      255 + (h[1]-255)*v,
      255 + (h[2]-255)*v,
      255 + (h[3]-255)*v
    )

    out[[length(out)+1]] <- as.data.table(shade)
    out[[length(out)+1]] <- as.data.table(tint)
  }
  dt <- rbindlist(out)
  setnames(dt, c("R","G","B"))
  unique_rgb_dt(dt)
}

make_primaries_secondaries <- function() {
  data.table(
    R=c(0,255,0,0,255,255,0,255),
    G=c(0,0,255,0,0,255,255,255),
    B=c(0,0,0,255,255,0,255,255)
  ) |> unique_rgb_dt()
}

# Sample points on the RGB cube faces (outer shell). Uses linear-light steps -> sRGB codes.
make_shell_faces <- function(k_face = 15L) {
  v <- linear_steps_srgb(as.integer(k_face))
  # 6 faces: R=0,R=255,G=0,G=255,B=0,B=255
  faces <- list(
    data.table(R=0L,   G=rep(v, each=length(v)), B=rep(v, times=length(v))),
    data.table(R=255L, G=rep(v, each=length(v)), B=rep(v, times=length(v))),
    data.table(G=0L,   R=rep(v, each=length(v)), B=rep(v, times=length(v))),
    data.table(G=255L, R=rep(v, each=length(v)), B=rep(v, times=length(v))),
    data.table(B=0L,   R=rep(v, each=length(v)), G=rep(v, times=length(v))),
    data.table(B=255L, R=rep(v, each=length(v)), G=rep(v, times=length(v)))
  )
  dt <- rbindlist(faces)
  setnames(dt, c("R","G","B"))
  unique_rgb_dt(dt)
}

# Off-grey rings (density-based): for each neutral v, for each ring magnitude, generate +/- axis offsets
# total_colors: FINAL chromatic target (after siblings) -> used to estimate sampling density
make_offgrey_rings <- function(neutrals, rings = 2L, total_colors = 6000L) {
  rings <- as.integer(rings)
  if (rings <= 0) return(data.table(R=integer(), G=integer(), B=integer()))
  nd <- nrow(neutrals)
  if (nd == 0) return(data.table(R=integer(), G=integer(), B=integer()))
  
  N <- max(1L, as.integer(total_colors))
  
  # Expected spacing scale in RGB codes (heuristic): ~255 / N^(1/3)
  base <- round(255 / (N^(1/3)))
  
  # Density-based max radius from expected RGB spacing, clamped to keep 3-ring output meaningful.
  delta_max <- max(3L, min(30L, as.integer(round(2 * base))))
  
  round_half_up <- function(x) floor(x + 0.5)

  # Start with evenly spaced radii
  mags <- as.integer(round_half_up(seq(1, delta_max, length.out = rings)))
  mags <- sort(unique(pmax(1L, mags)))

  # Force include endpoints (helps for small ring counts)
  mags <- sort(unique(c(1L, mags, delta_max)))

  # If we still have fewer than 'rings' distinct values, fill with missing ints
  if (length(mags) < rings) {
    missing <- setdiff(seq_len(delta_max), mags)
    mags <- sort(c(mags, head(missing, rings - length(mags))))
  }

  # If we have more than 'rings' (can happen after forcing endpoints), thin evenly
  if (length(mags) > rings) {
   keep_idx <- unique(as.integer(round_half_up(seq(1, length(mags), length.out = rings))))
   mags <- mags[keep_idx]
}
  
dirs <- rbind(
    c( 1, 0, 0), c(0,  1, 0), c(0, 0,  1),
    c(-1, 0, 0), c(0, -1, 0), c(0, 0, -1)
  )
  
  out <- vector("list", nd * length(mags) * nrow(dirs))
  idx <- 1L
  
  for (i in seq_len(nd)) {
    base_rgb <- as.integer(neutrals[i, .(R, G, B)])
    for (m in mags) {
      for (d in seq_len(nrow(dirs))) {
        rgb <- base_rgb + as.integer(dirs[d,] * m)
        out[[idx]] <- data.table(R=clip255(rgb[1]), G=clip255(rgb[2]), B=clip255(rgb[3]))
        idx <- idx + 1L
      }
    }
  }
  
  result <- unique(rbindlist(out))
  attr(result, "delta_max_used") <- delta_max
  result
}

# -------------------- fill: greedy maximin in OKLab --------------------
# Deterministic given seed.
# Strategy:
# 1) Build candidate RGBs via stratified jittered grid + random.
# 2) Convert to OKLab.
# 3) Greedy selection maximizing distance to nearest selected (maximin).
make_fill_poisson_oklab <- function(n_need, seed=12345, candidate_mult=6L, max_candidates=40000L) {
  n_need <- as.integer(n_need)
  if (n_need <= 0) return(data.table(R=integer(),G=integer(),B=integer()))
  set.seed(as.integer(seed))

  # Candidate count
  n_cand <- min(as.integer(max_candidates), as.integer(n_need * candidate_mult))
  n_cand <- max(n_cand, min(5000L, max_candidates))

  # Build candidates: mix of jittered grid and random
  # jittered grid
  m <- ceiling(n_cand^(1/3))
  grid_vals <- seq(0, 255, length.out = m)
  # sample m^3 points without building full cube
  idxs <- sample.int(m^3, size = min(n_cand, m^3), replace = FALSE)
  # map linear index -> (i,j,k)
  i <- ((idxs - 1) %% m) + 1
  j <- (((idxs - 1) %/% m) %% m) + 1
  k <- ((idxs - 1) %/% (m*m)) + 1
  cand <- cbind(grid_vals[i], grid_vals[j], grid_vals[k])

  # add jitter
  jitter <- matrix(runif(nrow(cand)*3, -0.45, 0.45), ncol=3)
  step <- if (m > 1) 255/(m-1) else 255
  cand <- cand + jitter * step
  cand <- clip255(cand)

  # ensure uniqueness
  cand_dt <- unique_rgb_dt(as.data.table(cand))
  cand <- as.matrix(cand_dt[, .(R,G,B)])
  n_cand <- nrow(cand)
  if (n_cand < n_need) {
    # fallback: add pure random candidates
    extra <- matrix(sample.int(256, size=(n_need - n_cand)*3, replace=TRUE)-1L, ncol=3)
    cand <- rbind(cand, extra)
    cand_dt <- unique_rgb_dt(as.data.table(cand))
    cand <- as.matrix(cand_dt[, .(R,G,B)])
    n_cand <- nrow(cand)
  }

  ok <- rgb255_to_oklab(cand)

  # Greedy maximin (maximise distance to nearest selected), without re-selecting the same candidate
  n_need <- min(n_need, n_cand)
  sel_idx <- integer(n_need)
  sel <- rep(FALSE, n_cand)

  sel_idx[1] <- sample.int(n_cand, 1)
  sel[sel_idx[1]] <- TRUE

  min_d2 <- rowSums((ok - ok[sel_idx[1],])^2)
  min_d2[sel] <- -Inf

  if (n_need >= 2) {
    for (s in 2:n_need) {
      jmax <- which.max(min_d2)
      if (!is.finite(min_d2[jmax])) break

      sel_idx[s] <- jmax
      sel[jmax] <- TRUE

      d2 <- rowSums((ok - ok[jmax,])^2)
      min_d2 <- pmin(min_d2, d2)
      min_d2[sel] <- -Inf
    }
  }

  sel_idx <- sel_idx[sel_idx > 0]
  fill <- cand[sel_idx, , drop=FALSE]
  unique_rgb_dt(as.data.table(fill))
}

make_fill_poisson_oklab_balanced <- function(n_need,
                                             seed = 12345,
                                             n_bins = 10L,
                                             bin_weights = c(2,4,7,10,13,15,15,13,11,11),  # <-- your tweak
                                             candidate_mult = 10L,
                                             max_candidates = 80000L) {
  n_need <- as.integer(n_need)
  if (n_need <= 0) return(data.table(R=integer(), G=integer(), B=integer()))
  set.seed(as.integer(seed))
  
  # --- Generate candidates in LINEAR RGB (fixes "too many brights") ---
  n_cand <- min(as.integer(max_candidates), as.integer(n_need * candidate_mult))
  n_cand <- max(n_cand, 15000L)
  
  linear_to_srgb01 <- function(u) {
    ifelse(u <= 0.0031308, 12.92*u, 1.055*(u^(1/2.4)) - 0.055)
  }
  
  lin <- matrix(runif(n_cand * 3L), ncol = 3L)
  srgb01 <- linear_to_srgb01(lin)
  cand <- clip255(round(srgb01 * 255))
  
  cand_dt <- unique_rgb_dt(as.data.table(cand))
  cand <- as.matrix(cand_dt[, .(R,G,B)])
  
  ok <- rgb255_to_oklab(cand)  # your existing converter
  L <- ok[,1]                  # OKLab L ~ 0..1
  
  # --- Bin by lightness, allocate picks per bin ---
  n_bins <- as.integer(n_bins)
  if (length(bin_weights) != n_bins) stop("bin_weights length must equal n_bins")
  
  w <- bin_weights / sum(bin_weights)
  
  bin_id <- pmin(n_bins, pmax(1L, as.integer(floor(L * n_bins)) + 1L))
  
  n_per <- pmax(0L, as.integer(round(w * n_need)))
  
  # fix rounding drift to hit n_need exactly
  diff <- n_need - sum(n_per)
  if (diff != 0) {
    order_bins <- order(abs((seq_len(n_bins) - 0.5) - n_bins/2))
    k <- 1L
    while (diff != 0) {
      bi <- order_bins[k]
      if (diff > 0) {
        n_per[bi] <- n_per[bi] + 1L
        diff <- diff - 1L
      } else if (n_per[bi] > 0) {
        n_per[bi] <- n_per[bi] - 1L
        diff <- diff + 1L
      }
      k <- if (k == length(order_bins)) 1L else (k + 1L)
    }
  }
  
  # --- Greedy maximin selection inside a candidate subset ---
  select_maximin <- function(idx, k) {
    if (!length(idx) || k <= 0) return(integer())
    k <- min(k, length(idx))
    
    sub_ok <- ok[idx, , drop=FALSE]
    n <- nrow(sub_ok)
    
    sel <- rep(FALSE, n)
    sel_i <- integer(k)
    
    sel_i[1] <- sample.int(n, 1)
    sel[sel_i[1]] <- TRUE
    
    min_d2 <- rowSums((sub_ok - sub_ok[sel_i[1],])^2)
    min_d2[sel] <- -Inf
    
    if (k >= 2) {
      for (s in 2:k) {
        j <- which.max(min_d2)
        if (!is.finite(min_d2[j])) break
        sel_i[s] <- j
        sel[j] <- TRUE
        
        d2 <- rowSums((sub_ok - sub_ok[j,])^2)
        min_d2 <- pmin(min_d2, d2)
        min_d2[sel] <- -Inf
      }
    }
    
    idx[sel_i[sel_i > 0]]
  }
  
  chosen <- integer()
  
  for (bi in seq_len(n_bins)) {
    idx <- which(bin_id == bi)
    chosen <- c(chosen, select_maximin(idx, n_per[bi]))
  }
  
  # Top-up if bins were too small (rare)
  if (length(chosen) < n_need) {
    remaining <- setdiff(seq_len(nrow(cand)), chosen)
    chosen <- c(chosen, select_maximin(remaining, n_need - length(chosen)))
  }
  
  fill <- cand[chosen, , drop=FALSE]
  unique_rgb_dt(as.data.table(fill))
}

# -------------------- layout / paging --------------------

get_page_mm <- function(paper, page_w_mm, page_h_mm) {
  if (paper == "A3")   return(list(w=297L, h=420L))
  if (paper == "A3+")  return(list(w=329L, h=483L))  # common 13x19"
  list(w=as.integer(page_w_mm), h=as.integer(page_h_mm))
}

compute_grid_paged <- function(page_w_mm, page_h_mm,
                               patch_mm, gap_mm,
                               margin_side_mm, margin_leading_mm, margin_trailing_mm,
                               max_chart_w_mm = 300,
                               n_patches,
                               max_pages = 3L) {
  usable_w_mm <- min(max_chart_w_mm, page_w_mm - 2 * margin_side_mm)
  usable_h_mm <- page_h_mm - margin_leading_mm - margin_trailing_mm
  pitch_mm <- patch_mm + gap_mm

  cols <- max(1L, floor((usable_w_mm + gap_mm) / pitch_mm))
  max_rows_fit <- max(1L, floor((usable_h_mm + gap_mm) / pitch_mm))

  total_rows_needed <- as.integer(ceiling(n_patches / cols))
  pages_needed <- as.integer(ceiling(total_rows_needed / max_rows_fit))
  pages <- max(1L, min(as.integer(max_pages), pages_needed))

  if (pages_needed > max_pages) {
    # won't fit without more pages
    return(list(fits=FALSE, cols=cols, max_rows_fit=max_rows_fit, pages_needed=pages_needed))
  }

  rows_per_page <- as.integer(ceiling(total_rows_needed / pages))
  rows_per_page <- min(rows_per_page, max_rows_fit)

  total_capacity <- cols * rows_per_page * pages

  list(
    fits=TRUE,
    cols=cols,
    rows_per_page=rows_per_page,
    pages=pages,
    per_page=cols*rows_per_page,
    total_capacity=total_capacity,
    max_rows_fit=max_rows_fit
  )
}

pad_to_capacity <- function(dt, capacity, mode=c("alternate","white","black")) {
  mode <- match.arg(mode)
  n <- nrow(dt)
  if (n >= capacity) return(dt)
  need <- capacity - n
  pad <- switch(
    mode,
    "white" = data.table(R=rep(255L, need), G=rep(255L, need), B=rep(255L, need)),
    "black" = data.table(R=rep(0L, need),   G=rep(0L, need),   B=rep(0L, need)),
    "alternate" = {
      s <- seq_len(need)
      w <- (s %% 2) == 1
      data.table(R=ifelse(w,255L,0L), G=ifelse(w,255L,0L), B=ifelse(w,255L,0L))
    }
  )
  rbindlist(list(dt, pad), use.names = TRUE)
}

add_white_black_repeats <- function(dt, n_each) {
  n_each <- as.integer(n_each)
  if (n_each <= 0) return(dt)
  wb <- rbindlist(list(
    data.table(R=rep(255L, n_each), G=rep(255L, n_each), B=rep(255L, n_each)),
    data.table(R=rep(0L,   n_each), G=rep(0L,   n_each), B=rep(0L,   n_each))
  ))
  rbindlist(list(dt, wb), use.names = TRUE)
}

# -------------------- outputs: CGATS + TIFF pages + MYIRO XML --------------------

write_cgats_rgb <- function(dt, path, title = "RGB_Target_PrintOrder") {
  n <- nrow(dt)
  lines <- c(
    "CGATS.17",
    paste0("ORIGINATOR\tRGB_Generator"),
    paste0("DESCRIPTOR\t", title),
    "NUMBER_OF_FIELDS\t4",
    paste0("NUMBER_OF_SETS\t", n),
    "BEGIN_DATA_FORMAT",
    "SAMPLE_ID\tRGB_R\tRGB_G\tRGB_B",
    "END_DATA_FORMAT",
    "BEGIN_DATA"
  )
  data_lines <- sprintf("%d\t%d\t%d\t%d", seq_len(n), dt$R, dt$G, dt$B)
  writeLines(c(lines, data_lines, "END_DATA"), con = path, useBytes = TRUE)
}

# MYIROtools chart definition XML modeled after sample:
# - pages via <Sheet Page="1"> ... </Sheet>
# - patches as <Patch Row=".." Column=".." ColorValue1="R*100" ...>
# Geometry in 0.01mm (integers). PatchColor values in 0..25500.
write_myiro_chart_xml <- function(dt, path,
                                  page_w_mm, page_h_mm,
                                  cols, rows_per_page, pages,
                                  patch_mm, gap_mm,
                                  margin_left_mm, margin_top_mm,
                                  title="RGB Target") {
  # MYIROtools ChartDefinition schema (v5), compatible with charts imported/created in MYIROtools.
  # Assumes dt is already in PRINT ORDER and padded so nrow(dt) == cols*rows_per_page*pages
  stopifnot(nrow(dt) == cols * rows_per_page * pages)
  
  mm100 <- function(x) as.integer(round(x * 100))
  
  # MYIROtools uses "Gap" fields as the PITCH (patch + gap), not the gap itself
  pitch <- mm100(patch_mm + gap_mm)      # PatchWidthGap / PatchHeightGap
  pw <- mm100(patch_mm)                  # PatchWidth
  ph <- mm100(patch_mm)                  # PatchHeight
  left <- mm100(margin_left_mm)          # LeftTopX
  top  <- mm100(margin_top_mm)           # LeftTopY
  
  # Sheet size should be the actual page size in 0.01mm
  sheetW <- mm100(page_w_mm)
  sheetH <- mm100(page_h_mm)
  
  total <- nrow(dt)
  per_page <- cols * rows_per_page
  
  # UUID-like token (MYIROtools often has one; format not strict)
  set.seed(1)
  uuid <- paste0(sample(c(letters[1:6], 0:9), 32, TRUE), collapse="")
  
  lines <- character()
  lines <- c(lines, '<?xml version="1.0" encoding="UTF-8"?>')
  lines <- c(lines, '<Data Identifier="ChartDefinition" Version="05.00" Application="MYIROtools">')
  
  # Keep the Chart attribute set broad and MYIRO-like (empty strings are acceptable)
  lines <- c(lines, sprintf(
    '  <Chart UNIT="0" PreviewProfile="" fdxMode="" fdxMode2="" CreateTime="" Kind="Normal" fd9RecogVerify="" Barcode="0" fd9Aperture="" i1proMode="" fd7Mode="" UUID="%s" ChartName="%s" PrintFileThumbnail="" Recog="Combine" PrintFile="" ChartPatchCount="%d" SpotScanLock="" MeasurementObjectType="" SheetCount="%d" i1pro2Mode="">',
    uuid, title, total, pages
  ))
  
  for (p in seq_len(pages)) {
    start <- (p-1)*per_page + 1
    end   <- p*per_page
    sub   <- dt[start:end]
    
    lines <- c(lines, sprintf(
      '    <Sheet AreaCount="1" BarcodeHeight="0" PatchCount="%d" Height="%d" BarcodeWidth="0" Page="%d" BarcodeLeftTopX="0" BarcodeLeftTopY="0" Width="%d">',
      per_page, sheetH, p, sheetW
    ))
    
    lines <- c(lines, sprintf(
      '      <Area PatchWidthGap="%d" PatchHeightGap="%d" PatchCount="%d" RowCount="%d" ColumnCount="%d" Number="1" PatchWidth="%d" PatchHeight="%d" LeftTopX="%d" LeftTopY="%d">',
      pitch, pitch, per_page, rows_per_page, cols, pw, ph, left, top
    ))
    
    # Sequential IDs (1..total) across pages
    base_id <- (p-1)*per_page
    idx <- 1L
    for (r in seq_len(rows_per_page)) {
      for (c in seq_len(cols)) {
        rgb <- sub[idx,]
        id  <- base_id + idx
        lines <- c(lines, sprintf(
          '        <Patch ColorValue1="%d" ID="%d" ColorValue2="%d" ColorValue3="%d" Attributes="" ColorFormat="RGB" Group="" Column="%d" Row="%d" />',
          as.integer(rgb$R*100), id, as.integer(rgb$G*100), as.integer(rgb$B*100), c, r
        ))
        idx <- idx + 1L
      }
    }
    
    lines <- c(lines, '      </Area>')
    lines <- c(lines, '    </Sheet>')
  }
  
  lines <- c(lines, '  </Chart>')
  lines <- c(lines, '</Data>')
  
  writeLines(lines, con = path, useBytes = TRUE)
}

write_tiff_pages_forced <- function(dt, out_dir,
                                    page_w_mm, page_h_mm,
                                    dpi,
                                    patch_mm, gap_mm,
                                    margin_left_mm, margin_top_mm, margin_trailing_mm,
                                    cols, rows_per_page, pages,
                                    filename_prefix="pages_") {
  mm_to_px <- function(mm) as.integer(round(mm * dpi / 25.4))
  page_w_px <- mm_to_px(page_w_mm)
  page_h_px <- mm_to_px(page_h_mm)
  patch_px  <- mm_to_px(patch_mm)
  gap_px    <- mm_to_px(gap_mm)
  left_px   <- mm_to_px(margin_left_mm)
  top_px    <- mm_to_px(margin_top_mm)
  trailing_px <- mm_to_px(margin_trailing_mm)
  y_limit <- page_h_px - trailing_px

  per_page <- cols * rows_per_page

  files <- character(pages)
  for (p in seq_len(pages)) {
    start <- (p-1)*per_page + 1
    end <- p*per_page
    sub <- dt[start:end]

    img <- array(255L, dim=c(page_h_px, page_w_px, 3))
    idx <- 1
    for (r in seq_len(rows_per_page)) {
      for (c in seq_len(cols)) {
        x <- left_px + (c-1)*(patch_px + gap_px) + 1
        y <- top_px  + (r-1)*(patch_px + gap_px) + 1
        x1 <- min(page_w_px, x + patch_px - 1)
        y1 <- min(y_limit,  y + patch_px - 1)

        img[y:y1, x:x1, 1] <- as.integer(sub$R[idx])
        img[y:y1, x:x1, 2] <- as.integer(sub$G[idx])
        img[y:y1, x:x1, 3] <- as.integer(sub$B[idx])
        idx <- idx + 1
      }
    }

    fname <- sprintf("%s%03d.tif", filename_prefix, p)
    fpath <- file.path(out_dir, fname)
    writeTIFF(img/255, fpath, compression="none")
    files[p] <- fpath
  }
  files
}

# -------------------- app --------------------

ui <- fluidPage(
  titlePanel("High-quality RGB Target Generator (MYIRO XML)"),
  sidebarLayout(
    sidebarPanel(
      numericInput("total_colors", "Chromatic target colors (final after siblings; before WB/pad)", value = 6000, min = 4000, max = 12000, step = 100),
      numericInput("neutral_steps", "Neutral greys count (R=G=B)", value = 33, min = 5, max = 256, step = 1),
      numericInput("offgrey_rings", "Off-grey ring count around neutral axis", value = 2, min = 0, max = 8, step = 1),

      checkboxInput("use_hue_ramps", "Include primary/secondary ramps (tint+shade)", value = TRUE),
      checkboxInput("add_siblings", "Add siblings (RGB cyclic permutations)", value = TRUE),

      tags$hr(),
      numericInput("wb_repeats_each", "Extra white patches AND extra black patches (each)", value = 24, min = 0, max = 5000, step = 1),

      tags$hr(),
      numericInput("seed", "Seed (reproducible)", value = 12345, min = 1, max = 9999999, step = 1),
      selectInput("pad_mode", "Padding color", choices = c("alternate"="alternate","white"="white","black"="black"), selected = "alternate"),

      tags$hr(),
      numericInput("max_pages", "Max pages allowed", value = 3, min = 1, max = 50, step = 1),
      numericInput("dpi", "TIFF DPI", value = 300, min = 150, max = 600, step = 10),
      selectInput("paper", "Paper", choices = c("A3"="A3","A3+"="A3+","Custom"="custom"), selected="A3+"),
      conditionalPanel(
        condition = "input.paper == 'custom'",
        numericInput("page_w_mm", "Page width (mm)", value = 329, min = 100, max = 2000, step = 1),
        numericInput("page_h_mm", "Page height (mm)", value = 483, min = 100, max = 2000, step = 1)
      ),
      numericInput("patch_mm", "Patch size (mm)", value = 6, min = 6, max = 30, step = 1),
      numericInput("gap_mm", "Gap (mm)", value = 1, min = 0, max = 5, step = 1),

      tags$hr(),
      downloadButton("download_zip", "Download ZIP (MYIRO XML + TIFF pages + CGATS + CSV)")
    ),
    mainPanel(
      verbatimTextOutput("stats"),
      tableOutput("preview")
    )
  )
)

server <- function(input, output, session) {

  # Core generation: returns a list with dt (print order, padded), gridp, and base stats
  generated <- reactive({

    set.seed(as.integer(input$seed))

    k <- as.integer(input$neutral_steps)
    rings <- as.integer(input$offgrey_rings)

    neutrals <- make_neutrals(k)
    primsec <- make_primaries_secondaries()
    ramps <- data.table(R=integer(),G=integer(),B=integer())
    if (isTRUE(input$use_hue_ramps)) ramps <- make_hue_ramps(k)

    offgrey <- make_offgrey_rings(
      neutrals,
      rings = input$offgrey_rings,
      total_colors = input$total_colors
    )
    delta_max_used <- attr(offgrey, "delta_max_used")

    mandatory <- unique(rbindlist(list(neutrals, primsec, ramps, offgrey)))

    # Apply siblings if enabled (mandatory too)
    if (isTRUE(input$add_siblings)) mandatory <- add_siblings(mandatory)

    # Build fill pool WITHOUT truncating mandatory
    target_chromatic <- as.integer(input$total_colors)

    # If mandatory already exceeds target, we KEEP it (quality > arbitrary count),
    # but we will still enforce max_pages via paging.
base <- mandatory

# --- Shell budget (fixed percentage) ---
shell_pct <- 0.20  # 20% of target goes to shell
target_chromatic <- as.integer(input$total_colors)
target_shell <- as.integer(round(shell_pct * target_chromatic))

# Generate shell candidates (faces) and take what we can use
shell_all <- make_shell_faces(k_face = 15L)

if (isTRUE(input$add_siblings)) shell_all <- add_siblings(shell_all)
shell_all <- shell_all[!base, on=.(R,G,B)]

# If shell_all is larger than needed, downsample deterministically
if (nrow(shell_all) > target_shell) {
  set.seed(as.integer(input$seed) + 77)
  shell_all <- shell_all[sample.int(nrow(shell_all), target_shell)]
}

base <- unique(rbindlist(list(base, shell_all)))
n_base <- nrow(base)

# Remaining budget becomes volumetric fill
n_need_final <- max(0L, target_chromatic - n_base)

fill <- data.table(R=integer(),G=integer(),B=integer())
    if (n_need_final > 0) {
      # generate only as many "root" colors as needed; siblings will be handled implicitly by root->triplet
      n_root <- if (isTRUE(input$add_siblings)) ceiling(n_need_final / 3) else n_need_final
      
      root <- make_fill_poisson_oklab_balanced(n_root, seed = as.integer(input$seed) + 101, candidate_mult=6L, max_candidates=40000L)
      
      if (isTRUE(input$add_siblings)) {
        fill <- add_siblings(root)
      } else {
        fill <- root
      }
      
      fill <- fill[!base, on=.(R,G,B)]
    }
    chromatic <- unique(rbindlist(list(base, fill)))
    
    # If still short, top-up similarly (same rule)
    tries <- 0
    while (nrow(chromatic) < target_chromatic && tries < 6) {
      tries <- tries + 1
      need2 <- target_chromatic - nrow(chromatic)
      n_root2 <- if (isTRUE(input$add_siblings)) ceiling(need2 / 3) else need2
      
      root2 <- make_fill_poisson_oklab_balanced(n_root2, seed = as.integer(input$seed) + 101 + tries, candidate_mult=7L, max_candidates=50000L)
      extra <- if (isTRUE(input$add_siblings)) add_siblings(root2) else root2
      extra <- extra[!chromatic, on=.(R,G,B)]
      chromatic <- unique(rbindlist(list(chromatic, extra)))
    }

    # If we still undershot because of dedupe/siblings collisions, top-up iteratively
    tries <- 0
    while (nrow(chromatic) < target_chromatic && tries < 6) {
      tries <- tries + 1
      extra_need <- target_chromatic - nrow(chromatic)
      extra <- make_fill_poisson_oklab_balanced(extra_need, seed = as.integer(input$seed) + 1101 + tries, candidate_mult=7L, max_candidates=50000L)
      if (isTRUE(input$add_siblings)) extra <- add_siblings(extra)
      extra <- extra[!chromatic, on=.(R,G,B)]
      chromatic <- unique(rbindlist(list(chromatic, extra)))
    }

    # We allow overshoot (keep it) as long as it fits in max_pages
    # Add repeated white/black measurements AFTER chromatic selection
    dt <- copy(chromatic)
    dt <- add_white_black_repeats(dt, input$wb_repeats_each)

    # Scramble print order to counter drift
    set.seed(as.integer(input$seed))
    dt <- dt[sample.int(nrow(dt))]

    # Paging / capacity
    pg <- get_page_mm(input$paper, input$page_w_mm, input$page_h_mm)

    # MYIRO scan margin handling: we reserve a safe left/top margin for grid placement.
    # (Keep your FD-S2w-like minimums; adjust if your MYIROtools template expects different.)
    margin_side_mm <- 4
    margin_top_mm <- 23
    margin_trailing_mm <- 33

    gridp <- compute_grid_paged(
      page_w_mm = pg$w,
      page_h_mm = pg$h,
      patch_mm = as.integer(input$patch_mm),
      gap_mm = as.integer(input$gap_mm),
      margin_side_mm = margin_side_mm,
      margin_leading_mm = margin_top_mm,
      margin_trailing_mm = margin_trailing_mm,
      max_chart_w_mm = 300,
      n_patches = nrow(dt),
      max_pages = as.integer(input$max_pages)
    )

    if (!isTRUE(gridp$fits)) {
      stop(sprintf("Does not fit within max_pages=%d. Need at least %d pages with current patch/gap and margins.",
                   as.integer(input$max_pages), as.integer(gridp$pages_needed)))
    }

    # If there is spare capacity within max_pages, keep overshoot; then pad to full capacity
    dt <- pad_to_capacity(dt, capacity = gridp$total_capacity, mode = input$pad_mode)

    # Return final
    list(
      dt = dt,
      gridp = gridp,
      pg = pg,
      delta_max = delta_max_used,
      mandatory_count = nrow(mandatory),
      chromatic_count = nrow(chromatic)
    )
  })

  output$stats <- renderPrint({
    g <- generated()
    dt <- g$dt
    gridp <- g$gridp
    cat("Chromatic set (after siblings, incl mandatory+fill):", g$chromatic_count, "\n")
    cat("Mandatory set size:", g$mandatory_count, "\n")
    cat("Printed patches (incl WB repeats + padding):", nrow(dt), "\n")
    cat("Off-grey delta_max used:", g$delta_max, "\n\n")

    cat("Page (mm):", g$pg$w, "x", g$pg$h, "\n")
    cat("Grid:", gridp$cols, "cols x", gridp$rows_per_page, "rows/page x", gridp$pages, "pages\n")
    cat("Capacity:", gridp$total_capacity, "patches\n\n")

    print(dt[, .(R_min=min(R), R_max=max(R), G_min=min(G), G_max=max(G), B_min=min(B), B_max=max(B))])
  })

  output$preview <- renderTable({
    head(generated()$dt, 30)
  })

  output$download_zip <- downloadHandler(
    filename = function() "rgb_target_myiro.zip",
    content = function(file) {
      g <- generated()
      dt <- g$dt
      gridp <- g$gridp
      pg <- g$pg

      tmpdir <- tempdir()
      outdir <- file.path(tmpdir, paste0("target_", as.integer(Sys.time())))
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

      # outputs
      csv_path   <- file.path(outdir, "colors_print_order.csv")
      cgats_path <- file.path(outdir, "target.cgats.txt")
      xml_path   <- file.path(outdir, "chart_definition.xml")

      fwrite(dt, csv_path)
      write_cgats_rgb(dt, cgats_path, title = "RGB_Target_PrintOrder")

      # MYIRO XML:
      write_myiro_chart_xml(
        dt = dt,
        path = xml_path,
        page_w_mm = pg$w,
        page_h_mm = pg$h,
        cols = gridp$cols,
        rows_per_page = gridp$rows_per_page,
        pages = gridp$pages,
        patch_mm = as.integer(input$patch_mm),
        gap_mm = as.integer(input$gap_mm),
        margin_left_mm = 4,
        margin_top_mm = 23,
        title = "RGB Target"
      )

      # TIFF pages (forced to same grid)
      tiff_files <- write_tiff_pages_forced(
        dt = dt, out_dir = outdir,
        page_w_mm = pg$w, page_h_mm = pg$h,
        dpi = as.integer(input$dpi),
        patch_mm = as.integer(input$patch_mm),
        gap_mm = as.integer(input$gap_mm),
        margin_left_mm = 4,
        margin_top_mm = 23,
        margin_trailing_mm = 33,
        cols = gridp$cols,
        rows_per_page = gridp$rows_per_page,
        pages = gridp$pages,
        filename_prefix = "pages_"
      )

      manifest <- file.path(outdir, "manifest.txt")
      writeLines(c(
        paste0("chromatic_count=", g$chromatic_count),
        paste0("mandatory_count=", g$mandatory_count),
        paste0("printed_patches_total=", nrow(dt)),
        paste0("page_w_mm=", pg$w),
        paste0("page_h_mm=", pg$h),
        paste0("cols=", gridp$cols),
        paste0("rows_per_page=", gridp$rows_per_page),
        paste0("pages=", gridp$pages),
        paste0("capacity=", gridp$total_capacity),
        paste0("seed=", input$seed),
        paste0("neutral_steps=", input$neutral_steps),
        paste0("offgrey_rings=", input$offgrey_rings),
        paste0("wb_repeats_each=", input$wb_repeats_each),
        paste0("pad_mode=", input$pad_mode)
      ), manifest, useBytes = TRUE)

      files_abs <- list.files(outdir, full.names = TRUE)
      zip::zipr(zipfile = file, files = files_abs)
    }
  )
}

if (!identical(Sys.getenv("RGB_TARGET_APP_DISABLE_AUTORUN"), "1")) shinyApp(ui, server)
