# RGB Calibration Target Generator
# - final total after siblings (before WB repeats + padding)
# - neutrals are pure R=G=B
# - off-gray count = neutrals * multiplier
# - ramps (tint+shade) with same steps as neutrals
# - scramble print order (seeded)
# - add extra white + extra black repeats
# - pad to full rows (myirotools friendly)
# - multi-page TIFF output + CGATS + CSV inside a ZIP
# - A3 / A3+ / Custom page size in mm

install.packages(c("shiny", "data.table", "tiff", "zip"))

library(shiny)
library(data.table)
library(tiff)
library(zip)

clip255 <- function(x) {
  d <- dim(x)                      # keep matrix/array dims if present
  y <- as.integer(round(x))
  y[y < 0L] <- 0L
  y[y > 255L] <- 255L
  if (!is.null(d)) dim(y) <- d     # restore dims
  y
}

unique_rgb <- function(mat) {
  dt <- as.data.table(mat)
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

even_steps <- function(k) {
  if (k <= 1) return(0L)
  clip255(round(seq(0, 255, length.out = k)))
}

make_neutrals <- function(k) {
  v <- even_steps(k)
  data.table(R=v, G=v, B=v)
}

# Off-grey: total_offgray points, balanced across +/- axes, magnitudes spread 1..delta_max
make_offgray <- function(neutrals_dt, total_offgray, delta_max) {
  if (total_offgray <= 0) return(data.table(R=integer(), G=integer(), B=integer()))
  
  dirs <- rbind(
    c( 1, 0, 0),
    c( 0, 1, 0),
    c( 0, 0, 1),
    c(-1, 0, 0),
    c( 0,-1, 0),
    c( 0, 0,-1)
  )
  nd <- nrow(neutrals_dt)
  
  mags <- unique(clip255(round(seq(1, delta_max, length.out = max(1, min(delta_max, 8))))))
  if (length(mags) == 0) mags <- 1L
  
  out <- vector("list", total_offgray)
  for (i in seq_len(total_offgray)) {
    base_idx <- ((i - 1) %% nd) + 1
    dir_idx  <- ((i - 1) %% nrow(dirs)) + 1
    mag_idx  <- ((i - 1) %% length(mags)) + 1
    
    base <- as.integer(neutrals_dt[base_idx, .(R,G,B)])
    d <- dirs[dir_idx, ]
    m <- mags[mag_idx]
    rgb <- base + as.integer(d * m)
    
    out[[i]] <- data.table(R=clip255(rgb[1]), G=clip255(rgb[2]), B=clip255(rgb[3]))
  }
  unique(rbindlist(out))
}

make_hue_ramps <- function(k) {
  hues <- list(
    R = c(255,   0,   0),
    G = c(  0, 255,   0),
    B = c(  0,   0, 255),
    C = c(  0, 255, 255),
    M = c(255,   0, 255),
    Y = c(255, 255,   0)
  )
  t <- seq(0, 1, length.out = k)
  
  out <- list()
  for (h in hues) {
    h <- as.numeric(h)
    tint <- cbind(
      255 + (h[1] - 255) * t,
      255 + (h[2] - 255) * t,
      255 + (h[3] - 255) * t
    )
    shade <- cbind(h[1] * t, h[2] * t, h[3] * t)
    out[[length(out)+1]] <- as.data.table(tint)
    out[[length(out)+1]] <- as.data.table(shade)
  }
  dt <- rbindlist(out)
  setnames(dt, c("R","G","B"))
  unique_rgb(as.matrix(dt))
}

# Halton fill for remainder
halton_1d <- function(n, base) {
  out <- numeric(n)
  for (i in seq_len(n)) {
    f <- 1; r <- 0; idx <- i
    while (idx > 0) {
      f <- f / base
      r <- r + f * (idx %% base)
      idx <- idx %/% base
    }
    out[i] <- r
  }
  out
}
halton_3d <- function(n) cbind(halton_1d(n,2), halton_1d(n,3), halton_1d(n,5))

make_uniform_fill <- function(n) {
  if (n <= 0) return(data.table(R=integer(), G=integer(), B=integer()))
  pts <- halton_3d(n)
  rgb <- clip255(round(pts * 255))
  unique_rgb(rgb)
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

pad_to_full_rows <- function(dt, cols, mode = c("alternate", "white", "black")) {
  mode <- match.arg(mode)
  n <- nrow(dt)
  rem <- n %% cols
  if (rem == 0) return(dt)
  need <- cols - rem
  
  pad <- switch(
    mode,
    "white" = data.table(R=rep(255L, need), G=rep(255L, need), B=rep(255L, need)),
    "black" = data.table(R=rep(0L, need),   G=rep(0L, need),   B=rep(0L, need)),
    "alternate" = {
      seqs <- seq_len(need)
      is_white <- seqs %% 2 == 1
      data.table(
        R=ifelse(is_white, 255L, 0L),
        G=ifelse(is_white, 255L, 0L),
        B=ifelse(is_white, 255L, 0L)
      )
    }
  )
  rbindlist(list(dt, pad), use.names = TRUE)
}

write_cgats_rgb <- function(dt, path, title = "RGB_Target_PrintOrder") {
  n <- nrow(dt)
  lines <- c(
    "CGATS.17",
    paste0("ORIGINATOR\t", "RGB_Generator"),
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

get_page_mm <- function(paper, page_w_mm, page_h_mm) {
  if (paper == "A3")   return(list(w=297L, h=420L))
  if (paper == "A3+")  return(list(w=329L, h=483L))  # common 13x19" stock
  list(w=as.integer(page_w_mm), h=as.integer(page_h_mm))
}

compute_grid <- function(page_w_mm, page_h_mm, patch_mm, gap_mm,
                         margin_side_mm, margin_leading_mm, margin_trailing_mm, max_chart_w_mm) {
  usable_w_mm <- min(max_chart_w_mm, page_w_mm - 2 * margin_side_mm)
  usable_h_mm <- page_h_mm - margin_leading_mm - margin_trailing_mm
  pitch_mm <- patch_mm + gap_mm
  cols <- max(1, floor((usable_w_mm + gap_mm) / pitch_mm))
  rows <- max(1, floor((usable_h_mm + gap_mm) / pitch_mm))
  list(cols=cols, rows=rows, per_page=cols*rows)
}

write_patch_pages_tiff_mm <- function(dt, out_dir,
                                      page_w_mm, page_h_mm,
                                      dpi = 300,
                                      patch_mm = 6,
                                      gap_mm = 1,
                                      margin_side_mm = 4,
                                      margin_leading_mm = 23,
                                      margin_trailing_mm = 33,
                                      max_chart_w_mm = 300,
                                      filename_prefix = "page_") {
  n <- nrow(dt)
  if (n == 0) return(list(files=character(), cols=0, rows=0, per_page=0, pages=0))
  
  grid <- compute_grid(page_w_mm, page_h_mm, patch_mm, gap_mm,
                       margin_side_mm, margin_leading_mm, margin_trailing_mm, max_chart_w_mm)
  cols <- grid$cols
  rows <- grid$rows
  per_page <- grid$per_page
  
  mm_to_px <- function(mm) as.integer(round(mm * dpi / 25.4))
  
  page_w_px <- mm_to_px(page_w_mm)
  page_h_px <- mm_to_px(page_h_mm)
  patch_px  <- mm_to_px(patch_mm)
  gap_px    <- mm_to_px(gap_mm)
  
  margin_side_px     <- mm_to_px(margin_side_mm)
  margin_leading_px  <- mm_to_px(margin_leading_mm)
  margin_trailing_px <- mm_to_px(margin_trailing_mm)
  
  x0 <- margin_side_px + 1
  y0 <- margin_leading_px + 1
  y_limit <- page_h_px - margin_trailing_px
  
  pages <- ceiling(n / per_page)
  files <- character(pages)
  
  for (p in seq_len(pages)) {
    start <- (p - 1) * per_page + 1
    end <- min(n, p * per_page)
    sub <- dt[start:end]
    
    img <- array(255L, dim = c(page_h_px, page_w_px, 3))
    
    for (i in seq_len(nrow(sub))) {
      r <- as.integer(sub$R[i]); g <- as.integer(sub$G[i]); b <- as.integer(sub$B[i])
      rr <- ((i - 1) %/% cols) + 1
      cc <- ((i - 1) %% cols) + 1
      
      x <- x0 + (cc - 1) * (patch_px + gap_px)
      y <- y0 + (rr - 1) * (patch_px + gap_px)
      
      x1 <- min(page_w_px, x + patch_px - 1)
      y1 <- min(y_limit, y + patch_px - 1)
      
      img[y:y1, x:x1, 1] <- r
      img[y:y1, x:x1, 2] <- g
      img[y:y1, x:x1, 3] <- b
    }
    
    fname <- sprintf("%s%03d.tif", filename_prefix, p)
    fpath <- file.path(out_dir, fname)
    writeTIFF(img / 255, fpath, compression = "none")
    files[p] <- fpath
  }
  
  list(files=files, cols=cols, rows=rows, per_page=per_page, pages=pages)
}

ui <- fluidPage(
  titlePanel("RGB Calibration Target Generator (MYIRO-friendly)"),
  sidebarLayout(
    sidebarPanel(
      numericInput("total_colors", "Total colors (FINAL after siblings; before WB/pad)", value = 6000, min = 4000, max = 12000, step = 100),
      numericInput("neutral_steps", "Neutral grays count (R=G=B)", value = 33, min = 5, max = 512, step = 1),
      
      checkboxInput("use_offgray", "Add off-gray set", value = TRUE),
      numericInput("offgray_multiplier", "Off-gray multiplier (x neutrals)", value = 1, min = 0, max = 20, step = 1),
      
      checkboxInput("use_hue_ramps", "Add primary/secondary ramps (tint+shade)", value = TRUE),
      checkboxInput("add_siblings", "Add siblings (RGB cyclic permutations)", value = TRUE),
      
      tags$hr(),
      numericInput("wb_repeats_each", "Extra white patches AND extra black patches (each)", value = 24, min = 0, max = 5000, step = 1),
      
      tags$hr(),
      numericInput("seed", "Shuffle seed (reproducible)", value = 12345, min = 1, max = 9999999, step = 1),
      
      selectInput("pad_mode", "Padding color for full rows",
                  choices = c("alternate white/black"="alternate", "white"="white", "black"="black"),
                  selected = "alternate"),
      
      tags$hr(),
      numericInput("dpi", "TIFF DPI", value = 300, min = 150, max = 600, step = 10),
      selectInput("paper", "Paper", choices = c("A3"="A3", "A3+"="A3+", "Custom"="custom"), selected = "A3+"),
      
      conditionalPanel(
        condition = "input.paper == 'custom'",
        numericInput("page_w_mm", "Page width (mm)", value = 329, min = 100, max = 2000, step = 1),
        numericInput("page_h_mm", "Page height (mm)", value = 483, min = 100, max = 2000, step = 1)
      ),
      
      numericInput("patch_mm", "Patch size (mm) (>= 6)", value = 6, min = 6, max = 30, step = 1),
      numericInput("gap_mm", "Gap (mm)", value = 1, min = 0, max = 5, step = 1),
      
      tags$hr(),
      downloadButton("download_zip", "Download ZIP (pages + CGATS + CSV)")
    ),
    mainPanel(
      verbatimTextOutput("stats"),
      tableOutput("preview")
    )
  )
)

server <- function(input, output, session) {
  
  generated <- reactive({
    k <- as.integer(input$neutral_steps)
    neutrals <- make_neutrals(k)
    
    # Off-grey strength tied to neutral spacing; small & sensible for neutrality
    neutral_step <- 255 / (k - 1)
    delta_max <- as.integer(max(1, min(12, round(neutral_step / 3))))
    
    offgray <- data.table(R=integer(), G=integer(), B=integer())
    if (isTRUE(input$use_offgray) && input$offgray_multiplier > 0) {
      offgray <- make_offgray(neutrals,
                              total_offgray = k * as.integer(input$offgray_multiplier),
                              delta_max = delta_max)
    }
    
    ramps <- data.table(R=integer(), G=integer(), B=integer())
    if (isTRUE(input$use_hue_ramps)) ramps <- make_hue_ramps(k)
    
    base <- unique(rbindlist(list(neutrals, offgray, ramps)))
    
    target_final <- as.integer(input$total_colors)
    
    # Ensure enough before siblings, then siblings, then trim to target_final
    target_pre <- if (isTRUE(input$add_siblings)) ceiling(target_final / 3) else target_final
    overshoot  <- ceiling(target_pre * 1.35)
    
    fill_n <- max(0L, overshoot - nrow(base))
    fill <- make_uniform_fill(fill_n)
    
    dt <- unique(rbindlist(list(base, fill)))
    if (isTRUE(input$add_siblings)) dt <- add_siblings(dt)
    
    dt <- dt[order(R, G, B)]
    if (nrow(dt) > target_final) dt <- dt[1:target_final]
    
    # Add repeated measurements of paper white and max black
    dt <- add_white_black_repeats(dt, input$wb_repeats_each)
    
    # Scramble print order (reproducible)
    set.seed(as.integer(input$seed))
    dt <- dt[sample.int(nrow(dt))]
    
    # Compute cols for padding using page size and layout
    pg <- get_page_mm(input$paper, input$page_w_mm, input$page_h_mm)
    grid <- compute_grid(pg$w, pg$h,
                         patch_mm = as.integer(input$patch_mm),
                         gap_mm = as.integer(input$gap_mm),
                         margin_side_mm = 4,
                         margin_leading_mm = 23,
                         margin_trailing_mm = 33,
                         max_chart_w_mm = 300)
    
    # Pad to full rows (myirotools expects full rows)
    dt <- pad_to_full_rows(dt, cols = grid$cols, mode = input$pad_mode)
    
    dt
  })
  
  output$stats <- renderPrint({
    dt <- generated()
    k <- as.integer(input$neutral_steps)
    neutral_step <- 255 / (k - 1)
    delta_max <- as.integer(max(1, min(12, round(neutral_step / 3))))
    
    pg <- get_page_mm(input$paper, input$page_w_mm, input$page_h_mm)
    grid <- compute_grid(pg$w, pg$h,
                         patch_mm = as.integer(input$patch_mm),
                         gap_mm = as.integer(input$gap_mm),
                         margin_side_mm = 4,
                         margin_leading_mm = 23,
                         margin_trailing_mm = 33,
                         max_chart_w_mm = 300)
    
    cat("Base chromatic target (after siblings; before WB/pad):", input$total_colors, "\n")
    cat("Final printed patches (includes WB repeats + padding):", nrow(dt), "\n\n")
    
    cat("Neutrals:", k, "\n")
    cat("Off-gray:", if (input$use_offgray) paste0(k * input$offgray_multiplier, " (delta_max=", delta_max, ")") else "off", "\n")
    cat("Hue ramps:", if (input$use_hue_ramps) "on" else "off", "\n")
    cat("Siblings:", if (input$add_siblings) "on" else "off", "\n")
    cat("Extra whites each:", input$wb_repeats_each, " | extra blacks each:", input$wb_repeats_each, "\n")
    cat("Shuffle seed:", input$seed, "\n")
    cat("Padding mode:", input$pad_mode, "\n\n")
    
    cat("Page (mm):", pg$w, "x", pg$h, "\n")
    cat("Grid:", grid$cols, "cols x", grid$rows, "rows =", grid$per_page, "patches/page\n")
    cat("Pages:", ceiling(nrow(dt) / grid$per_page), "\n\n")
    
    print(dt[, .(R_min=min(R), R_max=max(R), G_min=min(G), G_max=max(G), B_min=min(B), B_max=max(B))])
  })
  
  output$preview <- renderTable({ head(generated(), 25) })
  
  output$download_zip <- downloadHandler(
    filename = function() "calibration_target.zip",
    content = function(file) {
      dt <- generated()
      
      pg <- get_page_mm(input$paper, input$page_w_mm, input$page_h_mm)
      
      tmpdir <- tempdir()
      outdir <- file.path(tmpdir, paste0("target_", as.integer(Sys.time())))
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
      
      # CSV + CGATS in PRINT ORDER (includes repeats & padding)
      csv_path <- file.path(outdir, "colors_print_order.csv")
      cgats_path <- file.path(outdir, "target.cgats.txt")
      fwrite(dt, csv_path)
      write_cgats_rgb(dt, cgats_path, title = "RGB_Target_PrintOrder")
      
      # TIFF pages
      info <- write_patch_pages_tiff_mm(
        dt, out_dir = outdir,
        page_w_mm = pg$w, page_h_mm = pg$h,
        dpi = as.integer(input$dpi),
        patch_mm = as.integer(input$patch_mm),
        gap_mm = as.integer(input$gap_mm),
        margin_side_mm = 4,
        margin_leading_mm = 23,
        margin_trailing_mm = 33,  # FIXED
        max_chart_w_mm = 300,
        filename_prefix = "page_"
      )
      
      manifest <- file.path(outdir, "manifest.txt")
      writeLines(c(
        paste0("base_chromatic_target_after_siblings=", input$total_colors),
        paste0("neutrals=", input$neutral_steps),
        paste0("offgray_multiplier=", input$offgray_multiplier),
        paste0("extra_white_each=", input$wb_repeats_each),
        paste0("extra_black_each=", input$wb_repeats_each),
        paste0("pad_mode=", input$pad_mode),
        paste0("shuffle_seed=", input$seed),
        paste0("paper=", input$paper),
        paste0("page_w_mm=", pg$w),
        paste0("page_h_mm=", pg$h),
        paste0("dpi=", input$dpi),
        paste0("patch_mm=", input$patch_mm),
        paste0("gap_mm=", input$gap_mm),
        paste0("cols_per_page=", info$cols),
        paste0("rows_per_page=", info$rows),
        paste0("patches_per_page=", info$per_page),
        paste0("pages=", info$pages),
        paste0("total_printed_patches=", nrow(dt))
      ), manifest, useBytes = TRUE)
      
      oldwd <- getwd()
      setwd(outdir)
      on.exit(setwd(oldwd), add = TRUE)
      
      files_to_zip <- list.files(outdir, full.names = FALSE)
      zip::zipr(zipfile = file, files = files_to_zip)
    }
  )
}

shinyApp(ui, server)
