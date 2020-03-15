library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")

pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/test_image2.pdf")
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(20-16, 32+16), ylim = c(72-4, 76+4), expand = FALSE) + scale_x_continuous(breaks = seq(-32,40,8)) + scale_y_continuous(breaks = seq(64,80,2)) + geom_segment(data = df, aes(x = x, y = y, xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.1,"cm"))) + xlab("") + ylab("")
dev.off()

pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/test_image.pdf")
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-20, 28), ylim = c(62, 76), expand = FALSE) + scale_x_continuous(breaks = seq(-24,32,8)) + scale_y_continuous(breaks = seq(64,74,2)) + geom_rect(xmin=2, xmax = 8, ymin = 66, ymax = 72, color = "blue", fill = NA, size = 1)
dev.off()
