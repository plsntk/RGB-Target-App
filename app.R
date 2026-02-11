# RGB Calibration Target Generator (MYIROtools XML output)
# - total colors is final AFTER siblings/permutations, BEFORE WB repeats and padding
# - neutral count is pure R=G=B
# - off-gray count = neutrals * multiplier
# - primary/secondary ramps use same step count as neutrals
# - scramble print order (seeded)
# - add extra white patches and extra black patches (each)
# - pad so that every page has the same number of rows (MYIROtools friendly)
# - outputs ZIP: colors_print_order.csv, target.cgats.txt, chart_definition.xml, pages_###.tif

install.packages(c("shiny", "data.table", "tiff", "zip"))

library(shiny)
library(data.table)
library(tiff)
library(zip)

# ---------------- helpers ----------------

clip255 <- function(x) {
  d <- dim(x)
  y <- as.integer(round(x))
  y[y < 0L] <- 0L
  y[y > 255L] <- 255L
  if (!is.null(d)) dim(y) <- d
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

# Off-grey: balanced across +/- axes; magnitudes spread 1..delta_max
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

pad_n <- function(n, mode = c("alternate","white","black")) {
  mode <- match.arg(mode)
  if (n <= 0) return(data.table(R=integer(), G=integer(), B=integer()))
  switch(
    mode,
    "white" = data.table(R=rep(255L, n), G=rep(255L, n), B=rep(255L, n)),
    "black" = data.table(R=rep(0L, n),   G=rep(0L, n),   B=rep(0L, n)),
    "alternate" = {
      seqs <- seq_len(n)
      is_white <- seqs %% 2 == 1
      data.table(
        R=ifelse(is_white, 255L, 0L),
        G=ifelse(is_white, 255L, 0L),
        B=ifelse(is_white, 255L, 0L)
      )
    }
  )
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

# Compute paging so every page has the same rows_per_page (<= physical max)
compute_grid_paged <- function(page_w_mm, page_h_mm,
                              patch_mm, gap_mm,
                              margin_side_mm, margin_top_mm, margin_bottom_mm,
                              max_chart_w_mm = 300,
                              n_patches) {
  usable_w_mm <- min(max_chart_w_mm, page_w_mm - 2 * margin_side_mm)
  usable_h_mm <- page_h_mm - margin_top_mm - margin_bottom_mm
  pitch_mm <- patch_mm + gap_mm

  cols <- max(1L, floor((usable_w_mm + gap_mm) / pitch_mm))
  max_rows_fit <- max(1L, floor((usable_h_mm + gap_mm) / pitch_mm))

  total_rows_needed <- as.integer(ceiling(n_patches / cols))
  pages <- as.integer(ceiling(total_rows_needed / max_rows_fit))
  if (pages < 1L) pages <- 1L

  rows_per_page <- as.integer(ceiling(total_rows_needed / pages))
  if (rows_per_page > max_rows_fit) rows_per_page <- max_rows_fit

  capacity <- cols * rows_per_page * pages

  list(cols=cols, rows_per_page=rows_per_page, pages=pages, capacity=capacity,
       pitch_mm=pitch_mm, max_rows_fit=max_rows_fit)
}

# MYIROtools chart definition XML writer (based on uploaded "Chart 5207 Patches.xml")
# Geometry units: 0.01 mm (mm * 100)
# Color values: 0..25500 (channel * 100)
write_myiro_chart_xml <- function(dt, path,
                                 chart_name,
                                 page_w_mm, page_h_mm,
                                 cols, rows_per_page, pages,
                                 patch_mm, gap_mm,
                                 margin_left_mm, margin_top_mm) {
  n <- nrow(dt)
  per_page <- cols * rows_per_page
  expected <- per_page * pages
  if (n != expected) stop("XML layout mismatch: patch count != cols*rows*pages")

  mm100 <- function(x) as.character(as.integer(round(x * 100)))
  ch100 <- function(x) as.character(as.integer(round(x * 100)))

  pitch_mm <- patch_mm + gap_mm

  # Header
  lines <- c(
    '<?xml version="1.0" encoding="utf-8"?>',
    '<Data Identifier="ChartDefinition" Version="05.00" Application="MYIROtools">',
    sprintf('<Chart UNIT="0" PreviewProfile="" fdxMode="" CreateTime="" Kind="Normal" fdxMode2="" fd9RecogVerify="" Barcode="0" fd9Aperture="" i1proMode="" fd7Mode="" UUID="" ChartName="%s" PrintFileThumbnail="" Recog="Combine" PrintFile="" ChartPatchCount="%d" SpotScanLock="" MeasurementObjectType="" SheetCount="%d" i1pro2Mode="">',
            gsub('"','&quot;', chart_name), n, pages)
  )

  # Sheets + Areas
  id_counter <- 1L
  for (p in seq_len(pages)) {
    # whole page sheet
    sheet_width  <- mm100(page_w_mm)
    sheet_height <- mm100(page_h_mm)

    lines <- c(lines,
      sprintf('<Sheet AreaCount="1" BarcodeHeight="0" PatchCount="%d" Height="%s" BarcodeWidth="0" Page="%d" BarcodeLeftTopX="0" BarcodeLeftTopY="0" Width="%s">',
              per_page, sheet_height, p, sheet_width),
      sprintf('<Area PatchWidthGap="%s" PatchHeightGap="%s" PatchCount="%d" RowCount="%d" ColumnCount="%d" Number="1" PatchWidth="%s" PatchHeight="%s" LeftTopX="%s" LeftTopY="%s">',
              mm100(pitch_mm), mm100(pitch_mm), per_page, rows_per_page, cols,
              mm100(patch_mm), mm100(patch_mm), mm100(margin_left_mm), mm100(margin_top_mm))
    )

    # patches for this page
    start <- (p - 1L) * per_page + 1L
    end   <- p * per_page
    sub <- dt[start:end]

    # order: row-major within area
    for (r in seq_len(rows_per_page)) {
      for (c in seq_len(cols)) {
        idx <- (r - 1L) * cols + c
        rgb <- sub[idx]
        lines <- c(lines,
          sprintf('<Patch ColorValue1="%s" ID="%d" ColorValue2="%s" ColorValue3="%s" Attributes="" ColorFormat="RGB" Group="" Column="%d" Row="%d" />',
                  ch100(rgb$R), id_counter, ch100(rgb$G), ch100(rgb$B), c, r)
        )
        id_counter <- id_counter + 1L
      }
    }

    lines <- c(lines, '</Area>', '</Sheet>')
  }

  lines <- c(lines, '</Chart>', '</Data>')
  writeLines(lines, path, useBytes = TRUE)
}

# Multi-page TIFF writer forced to cols/rows_per_page/pages
write_patch_pages_tiff_mm_fixed <- function(dt, out_dir,
                                           page_w_mm, page_h_mm,
                                           dpi,
                                           patch_mm, gap_mm,
                                           margin_left_mm, margin_top_mm, margin_bottom_mm,
                                           cols, rows_per_page, pages,
                                           filename_prefix = "pages_") {
  n <- nrow(dt)
  per_page <- cols * rows_per_page
  expected <- per_page * pages
  if (n != expected) stop("TIFF layout mismatch: patch count != cols*rows*pages")

  mm_to_px <- function(mm) as.integer(round(mm * dpi / 25.4))

  page_w_px <- mm_to_px(page_w_mm)
  page_h_px <- mm_to_px(page_h_mm)

  patch_px <- mm_to_px(patch_mm)
  gap_px   <- mm_to_px(gap_mm)
  pitch_px <- patch_px + gap_px

  x0 <- mm_to_px(margin_left_mm) + 1L
  y0 <- mm_to_px(margin_top_mm) + 1L
  y_limit <- page_h_px - mm_to_px(margin_bottom_mm)

  files <- character(pages)

  for (p in seq_len(pages)) {
    start <- (p - 1L) * per_page + 1L
    end   <- p * per_page
    sub <- dt[start:end]

    img <- array(255L, dim = c(page_h_px, page_w_px, 3))

    for (i in seq_len(per_page)) {
      rr <- ((i - 1L) %/% cols) + 1L
      cc <- ((i - 1L) %% cols) + 1L

      x <- x0 + (cc - 1L) * pitch_px
      y <- y0 + (rr - 1L) * pitch_px

      x1 <- min(page_w_px, x + patch_px - 1L)
      y1 <- min(y_limit,  y + patch_px - 1L)

      img[y:y1, x:x1, 1] <- sub$R[i]
      img[y:y1, x:x1, 2] <- sub$G[i]
      img[y:y1, x:x1, 3] <- sub$B[i]
    }

    fname <- sprintf("%s%03d.tif", filename_prefix, p)
    fpath <- file.path(out_dir, fname)
    writeTIFF(img / 255, fpath, compression = "none")
    files[p] <- fpath
  }

  files
}

# ---------------- UI ----------------

ui <- fluidPage(
  titlePanel("RGB Calibration Target Generator (MYIROtools XML)"),
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

      selectInput("pad_mode", "Padding color",
                  choices = c("alternate white/black"="alternate", "white"="white", "black"="black"),
                  selected = "alternate"),

      tags$hr(),
      textInput("chart_name", "Chart name", value = "RGB Target"),

      tags$hr(),
      numericInput("dpi", "TIFF DPI", value = 300, min = 150, max = 600, step = 10),
      selectInput("paper", "Paper", choices = c("A3"="A3", "A3+"="A3+", "Custom"="custom"), selected = "A3+"),
      conditionalPanel(
        condition = "input.paper == 'custom'",
        numericInput("page_w_mm", "Page width (mm)", value = 329, min = 100, max = 2000, step = 1),
        numericInput("page_h_mm", "Page height (mm)", value = 483, min = 100, max = 2000, step = 1)
      ),

      numericInput("patch_mm", "Patch size (mm)", value = 6, min = 2, max = 30, step = 1),
      numericInput("gap_mm", "Gap between patches (mm)", value = 0, min = 0, max = 10, step = 1),

      tags$hr(),
      numericInput("margin_left_mm", "Left/Right margin (mm)", value = 4, min = 0, max = 50, step = 1),
      numericInput("margin_top_mm", "Top margin (mm)", value = 23, min = 0, max = 80, step = 1),
      numericInput("margin_bottom_mm", "Bottom margin (mm)", value = 33, min = 0, max = 120, step = 1),

      tags$hr(),
      downloadButton("download_zip", "Download ZIP (TIFF + CGATS + MYIRO XML)")
    ),
    mainPanel(
      verbatimTextOutput("stats"),
      tableOutput("preview")
    )
  )
)

# ---------------- server ----------------

server <- function(input, output, session) {

  generated <- reactive({
    k <- as.integer(input$neutral_steps)
    neutrals <- make_neutrals(k)

    neutral_step <- 255 / (k - 1)
    delta_max <- as.integer(max(1, min(12, round(neutral_step / 3))))

    offgray <- data.table(R=integer(), G=integer(), B=integer())
    if (isTRUE(input$use_offgray) && input$offgray_multiplier > 0) {
      offgray <- make_offgray(neutrals, total_offgray = k * as.integer(input$offgray_multiplier), delta_max = delta_max)
    }

    ramps <- data.table(R=integer(), G=integer(), B=integer())
    if (isTRUE(input$use_hue_ramps)) ramps <- make_hue_ramps(k)

    base <- unique(rbindlist(list(neutrals, offgray, ramps)))

    target_final <- as.integer(input$total_colors)

    # Build enough pre-siblings, then siblings, then trim
    target_pre <- if (isTRUE(input$add_siblings)) ceiling(target_final / 3) else target_final
    overshoot  <- ceiling(target_pre * 1.35)

    fill_n <- max(0L, overshoot - nrow(base))
    fill <- make_uniform_fill(fill_n)

    dt <- unique(rbindlist(list(base, fill)))
    if (isTRUE(input$add_siblings)) dt <- add_siblings(dt)

    dt <- dt[order(R, G, B)]
    if (nrow(dt) > target_final) dt <- dt[1:target_final]

    # Add repeated paper white + black measurements
    dt <- add_white_black_repeats(dt, input$wb_repeats_each)

    # Scramble order (reproducible)
    set.seed(as.integer(input$seed))
    dt <- dt[sample.int(nrow(dt))]

    # Compute page/grid and pad so every page has equal rows
    pg <- get_page_mm(input$paper, input$page_w_mm, input$page_h_mm)
    gridp <- compute_grid_paged(
      page_w_mm = pg$w, page_h_mm = pg$h,
      patch_mm = as.integer(input$patch_mm),
      gap_mm   = as.integer(input$gap_mm),
      margin_side_mm = as.integer(input$margin_left_mm),
      margin_top_mm  = as.integer(input$margin_top_mm),
      margin_bottom_mm = as.integer(input$margin_bottom_mm),
      max_chart_w_mm = 300,
      n_patches = nrow(dt)
    )

    need <- gridp$capacity - nrow(dt)
    if (need > 0) dt <- rbindlist(list(dt, pad_n(need, mode = input$pad_mode)), use.names = TRUE)

    dt
  })

  output$stats <- renderPrint({
    dt <- generated()
    pg <- get_page_mm(input$paper, input$page_w_mm, input$page_h_mm)

    gridp <- compute_grid_paged(
      page_w_mm = pg$w, page_h_mm = pg$h,
      patch_mm = as.integer(input$patch_mm),
      gap_mm   = as.integer(input$gap_mm),
      margin_side_mm = as.integer(input$margin_left_mm),
      margin_top_mm  = as.integer(input$margin_top_mm),
      margin_bottom_mm = as.integer(input$margin_bottom_mm),
      max_chart_w_mm = 300,
      n_patches = nrow(dt)
    )

    cat("Base chromatic target (after siblings; before WB/pad):", input$total_colors, "\n")
    cat("Final printed patches (includes WB repeats + padding):", nrow(dt), "\n\n")
    cat("Page (mm):", pg$w, "x", pg$h, "\n")
    cat("Patch (mm):", input$patch_mm, " Gap (mm):", input$gap_mm, "\n")
    cat("Margins (mm): L/R", input$margin_left_mm, " Top", input$margin_top_mm, " Bottom", input$margin_bottom_mm, "\n\n")
    cat("Grid:", gridp$cols, "cols x", gridp$rows_per_page, "rows per page x", gridp$pages, "pages\n")
    cat("Capacity:", gridp$capacity, " patches\n")
    cat("Max rows that fit physically:", gridp$max_rows_fit, "\n\n")

    print(dt[, .(R_min=min(R), R_max=max(R), G_min=min(G), G_max=max(G), B_min=min(B), B_max=max(B))])
  })

  output$preview <- renderTable({ head(generated(), 25) })

  output$download_zip <- downloadHandler(
    filename = function() "calibration_target.zip",
    content = function(file) {
      dt <- generated()
      pg <- get_page_mm(input$paper, input$page_w_mm, input$page_h_mm)

      gridp <- compute_grid_paged(
        page_w_mm = pg$w, page_h_mm = pg$h,
        patch_mm = as.integer(input$patch_mm),
        gap_mm   = as.integer(input$gap_mm),
        margin_side_mm = as.integer(input$margin_left_mm),
        margin_top_mm  = as.integer(input$margin_top_mm),
        margin_bottom_mm = as.integer(input$margin_bottom_mm),
        max_chart_w_mm = 300,
        n_patches = nrow(dt)
      )

      # Safety: ensure dt exactly matches capacity
      need <- gridp$capacity - nrow(dt)
      if (need > 0) dt <- rbindlist(list(dt, pad_n(need, mode = input$pad_mode)), use.names = TRUE)
      if (nrow(dt) != gridp$capacity) stop("Internal error: dt != capacity")

      tmpdir <- tempdir()
      outdir <- file.path(tmpdir, paste0("target_", as.integer(Sys.time())))
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

      # CSV + CGATS in print order
      csv_path <- file.path(outdir, "colors_print_order.csv")
      cgats_path <- file.path(outdir, "target.cgats.txt")
      fwrite(dt, csv_path)
      write_cgats_rgb(dt, cgats_path, title = input$chart_name)

      # MYIROtools chart definition XML (page-aware)
      xml_path <- file.path(outdir, "chart_definition.xml")
      write_myiro_chart_xml(
        dt = dt, path = xml_path,
        chart_name = input$chart_name,
        page_w_mm = pg$w, page_h_mm = pg$h,
        cols = gridp$cols, rows_per_page = gridp$rows_per_page, pages = gridp$pages,
        patch_mm = as.integer(input$patch_mm), gap_mm = as.integer(input$gap_mm),
        margin_left_mm = as.integer(input$margin_left_mm),
        margin_top_mm  = as.integer(input$margin_top_mm)
      )

      # TIFF pages (match the same grid)
      write_patch_pages_tiff_mm_fixed(
        dt = dt, out_dir = outdir,
        page_w_mm = pg$w, page_h_mm = pg$h,
        dpi = as.integer(input$dpi),
        patch_mm = as.integer(input$patch_mm),
        gap_mm = as.integer(input$gap_mm),
        margin_left_mm = as.integer(input$margin_left_mm),
        margin_top_mm  = as.integer(input$margin_top_mm),
        margin_bottom_mm = as.integer(input$margin_bottom_mm),
        cols = gridp$cols, rows_per_page = gridp$rows_per_page, pages = gridp$pages,
        filename_prefix = "pages_"
      )

      # Manifest
      manifest <- file.path(outdir, "manifest.txt")
      writeLines(c(
        paste0("base_chromatic_target_after_siblings=", input$total_colors),
        paste0("final_printed_patches=", nrow(dt)),
        paste0("chart_name=", input$chart_name),
        paste0("page_w_mm=", pg$w),
        paste0("page_h_mm=", pg$h),
        paste0("dpi=", input$dpi),
        paste0("patch_mm=", input$patch_mm),
        paste0("gap_mm=", input$gap_mm),
        paste0("margin_left_mm=", input$margin_left_mm),
        paste0("margin_top_mm=", input$margin_top_mm),
        paste0("margin_bottom_mm=", input$margin_bottom_mm),
        paste0("cols=", gridp$cols),
        paste0("rows_per_page=", gridp$rows_per_page),
        paste0("pages=", gridp$pages),
        paste0("capacity=", gridp$capacity),
        paste0("seed=", input$seed)
      ), manifest, useBytes = TRUE)

      # Zip all files (absolute paths; no setwd)
      files_abs <- list.files(outdir, full.names = TRUE)
      zip::zipr(zipfile = file, files = files_abs)
    }
  )
}

shinyApp(ui, server)
