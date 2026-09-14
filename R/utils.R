#' Save a plot to PNG/PDF/TIFF
#'
#' @description
#' Save a ggplot (or cowplot) object to disk in PNG, PDF and TIFF under
#' `Fig/<FORMAT>/`. Mirrors the SILAC pipeline helper (supports `target_dir`).
#'
#' @param object Plot object to save.
#' @param filenames Base file name (without extension).
#' @param width Figure width in inches (required).
#' @param height Figure height in inches (required).
#' @param dpi Resolution for raster formats (default 600).
#' @param target_dir Optional subdirectory prepended to `Fig/`.
#'
#' @return Invisibly returns the file paths written.
#' @export
saveplot <- function(object,
                     filenames = filename,
                     width = NULL,
                     height = NULL,
                     dpi = 600,
                     target_dir = NULL) {
  
  if (is.null(width) | is.null(height)) {
    stop("Must input figure width and height")
  }
  
  if (is.null(filenames) || !nzchar(filenames)) {
    stop("Must input filenames")
  }
  
  if (is.null(target_dir) || !nzchar(target_dir)) {
    target_dir <- ''
  }
  
  # 主输出目录：Fig_{target_dir}
  fig_dir <- paste0("Fig_", target_dir)
  
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  purrr::map(c("TIFF", "PNG", "PDF"), function(x) {
    sub_dir <- file.path(fig_dir, x)
    if (!dir.exists(sub_dir)) {
      dir.create(sub_dir)
    }
  })
  
  pdf_filenames <- file.path(fig_dir, "PDF", paste0(filenames, ".pdf"))
  pdf_width <- width * 2
  pdf_height <- height * 2
  pdf_res <- dpi
  
  tiff_filenames <- file.path(fig_dir, "TIFF", paste0(filenames, ".tiff"))
  tiff_width <- width * 1000
  tiff_height <- height * 1000
  tiff_res <- dpi
  
  png_filenames <- file.path(fig_dir, "PNG", paste0(filenames, ".png"))
  png_width <- width * 500
  png_height <- height * 500
  
  if (dpi > 600) {
    png_res <- 300
  } else {
    png_res <- dpi / 2
  }
  
  if (class(object)[1] == "gList") {
    
    cairo_pdf(
      filename = pdf_filenames,
      width = pdf_width,
      height = pdf_height,
      fallback_resolution = pdf_res
    )
    grid::grid.draw(object)
    dev.off()
    
    tiff(
      filename = tiff_filenames,
      width = tiff_width,
      height = tiff_height,
      units = "px",
      res = tiff_res,
      compression = "lzw"
    )
    grid::grid.draw(object)
    dev.off()
    
    png(
      filename = png_filenames,
      width = png_width,
      height = png_height,
      res = png_res,
      units = "px"
    )
    grid::grid.draw(object)
    dev.off()
    
  } else {
    
    cairo_pdf(
      filename = pdf_filenames,
      width = pdf_width,
      height = pdf_height,
      fallback_resolution = pdf_res
    )
    print(object)
    dev.off()
    
    tiff(
      filename = tiff_filenames,
      width = tiff_width,
      height = tiff_height,
      units = "px",
      res = tiff_res,
      compression = "lzw"
    )
    print(object)
    dev.off()
    
    png(
      filename = png_filenames,
      width = png_width,
      height = png_height,
      res = png_res,
      units = "px"
    )
    print(object)
    dev.off()
  }
}

#' Remove saved plot files
#'
#' @description Remove PNG/PDF/TIFF files previously written by [saveplot()].
#'
#' @param filename Base file name (without extension).
#'
#' @return Invisibly NULL.
#' @export
remove_plot <- function(filename){
  if (file.exists(glue::glue('Fig/PDF/{filename}.pdf'))) {
    file.remove(glue::glue('Fig/PDF/{filename}.pdf'))
    cat('PDF removed\n')
  }

  if (file.exists(glue::glue('Fig/PNG/{filename}.png'))) {
    file.remove(glue::glue('Fig/PNG/{filename}.png'))
    cat('PNG removed\n')
  }
  if (file.exists(glue::glue('Fig/TIFF/{filename}.tiff'))) {
    file.remove(glue::glue('Fig/TIFF/{filename}.tiff'))
    cat('TIFF removed\n')
  }
}

#' Save a plot to PNG only
#'
#' @param object Plot object.
#' @param filenames Base file name.
#' @param width Width in inches (required).
#' @param height Height in inches (required).
#' @param dpi Resolution (default 600).
#'
#' @return Invisibly the written path.
#' @export
savepng <- function(object,filenames = filename,width = NULL,height = NULL,dpi = 600){
  if (is.null(width)| is.null(height)) {
    stop('Must input figure width and height')
  }
  
  if(!dir.exists("Fig")){
    dir.create("Fig")
  }
  map(c("TIFF","PNG","PDF"),function(x){
    if(!dir.exists(paste0('Fig/',x))){
      dir.create(paste0('Fig/',x))
    }
  })
  
  pdf_filenames = paste0('Fig/PDF/',filenames,'.pdf')
  pdf_width = width*2
  pdf_height = height*2
  pdf_res = dpi
  
  tiff_filenames = paste0('Fig/TIFF/',filenames,'.tiff')
  tiff_width = width*1000
  tiff_height = height*1000
  tiff_res = dpi
  
  png_filenames = paste0('Fig/PNG/',filenames,'.png')
  png_width = width*500
  png_height = height*500
  if (dpi > 600) {
    png_res = 300
  }else{
    png_res = dpi/2
  }
  if (class(object)[1] == 'gList') {
    
    
    png(filename = png_filenames,width = png_width, height = png_height, res = png_res,units = "px")
    grid.draw(object)
    dev.off()
  }else{
    
    
    png(filename = png_filenames,width = png_width, height = png_height, res = png_res,units = "px")
    print(object)
    dev.off()
  }
  
}

#' Print a timestamped message
#'
#' @param ... Passed to `cat()`.
#' @param verbose Logical; if FALSE the message is suppressed.
#'
#' @return Invisibly NULL.
#' @export
msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    message(glue::glue(...))
  }
}
