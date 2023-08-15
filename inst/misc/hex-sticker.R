## Author: Balasubramanian Narasimhan apres Abebe Geletu slides
## https://www.tu-ilmenau.de/fileadmin/Bereiche/IA/prozessoptimierung/vorlesungsskripte/abebe_geletu/IPM_Slides.pdf

library(hexSticker)
library(ggplot2)

gen_plot <- function() {
  p <- ggplot() + 
    theme_void() + 
    coord_fixed(ratio = 1) + 
    xlim(-2, 4) + 
    ylim(-2, 3)

  ## Add x and y axis
  p <- p +
    geom_segment(aes(x=-2, y=0, xend=3, yend=0), linewidth = 0.25, arrow=arrow(type="closed", angle = 10, length = unit(0.05, "inches"))) + 
    geom_segment(aes(x=0, y=-1.5, xend=0, yend=2.5), linewidth = 0.25, arrow=arrow(type="closed", angle = 10, length = unit(0.05, "inches")))

  rotate_point <- function(x, y, cx, cy, angle_deg) {
    angle_rad <- angle_deg * (pi / 180)
    x_new <- cx + (x - cx) * cos(angle_rad) - (y - cy) * sin(angle_rad)
    y_new <- cy + (x - cx) * sin(angle_rad) + (y - cy) * cos(angle_rad)
    return(c(x_new, y_new))
  }

  # Compute the rotated coordinates of the rectangle around the point (0.5, 0.5)
  bottom_left <- rotate_point(-1.0, -0.5, 0.5, 0.5, 30)
  top_left <- rotate_point(-1.0, 1.5, 0.5, 0.5, 30)
  top_right <- rotate_point(3.2, 1.5, 0.5, 0.5, 30)
  bottom_right <- rotate_point(3.2, -0.5, 0.5, 0.5, 30)

  rotated_rect <- data.frame(
    x=c(bottom_left[1], bottom_right[1], top_right[1], top_left[1]),
    y=c(bottom_left[2], bottom_right[2], top_right[2], top_left[2]) 
  )

  p <- p + geom_polygon(data=rotated_rect, aes(x, y), fill="darkcyan", alpha=0.8)

  # Label the rectangle
  p <- p + annotate("text", x=2.6, y=2, label="h(x) = 0", hjust="right", fontface = 'italic', size = unit(7, "pt"))

  # Add the g(x) region
  gregion <- data.frame(
    x=c(-0.5, 1, 1.5, 2.0, 2.2, 1),
    y=c(1, -1, -0.8, -0.1, 1.5, 1.8)
  )

  p <- p + geom_polygon(data=gregion, aes(x, y), fill="pink", alpha=0.65, colour="pink", linewidth=.75, linetype="solid")

  # Label g(x) region
  p <- p + annotate("text", x=1.75, y=-0.25, label="g(x) ≤ 0", hjust="right", fontface = 'italic', size = unit(7, "pt"))

  # Add the red lines
  redline <- data.frame(
    x=c(0, 0, 0.25, 1.633975, 2.0419426, 2.2, 1, 0),
    y=c(1.26667, 0.33333, 0, 0, 0.2355404, 1.5, 1.8, 1.26667)
  )

  p <- p + geom_path(data=redline, aes(x, y), color="red", linewidth=0.75)

  # Add the yellow path
  yellow_path <- data.frame(
    x=c(0.5, 0.8, 0.6, 1.0, 0.8, 1.2, 1.0, 1.3, 1.5),
    y=c(0.3, 0.6, 0.9, 1.2, 1.3, 1.4, 1.5, 1.6, 1.675)
  )

  p <- p + geom_path(data=yellow_path, aes(x, y), color="yellow", linewidth=0.5, lineend="round")

  # Add the circle at the end of the yellow path
  p <- p + geom_point(data=data.frame(x=1.5, y=1.675), aes(x, y), shape=21, size=2, color="black", fill="black")

  # Display the plot
  p
}

s <- sticker(gen_plot(),
             package="piqp", p_color = "black", p_size=20,
             p_y = 1.6,
             s_x=1.0, s_y=.8, s_width=2.0, s_height=1.8,
             h_fill = "white", h_color = "darkred",
             filename="logo.png")
