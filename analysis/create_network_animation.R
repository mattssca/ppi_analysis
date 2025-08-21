#load pacakges
library(magick)
library(pdftools)

#list all PDF files in the target directory (change path as needed)
pdf_files <- list.files(path = "viz/networks/test/", pattern = "*.pdf", full.names = TRUE)

#convert PDFs to PNGs
png_files <- character(length(pdf_files))
for (i in seq_along(pdf_files)) {
  img <- image_read_pdf(pdf_files[i], density = 300)
  png_file <- sub(".pdf$", ".png", pdf_files[i])
  image_write(img, path = png_file, format = "png")
  png_files[i] <- png_file
}

#read PNGs and create animation
img_list <- image_read(png_files)
animation <- image_animate(image_join(img_list), fps = 1) # 1 frame per second

#save animated GIF
image_write(animation, path = "viz/networks/network_animation.gif")
