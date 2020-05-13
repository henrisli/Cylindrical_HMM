#https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html

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

pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/map_of_area.pdf")
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-13, 22), ylim = c(62, 76), expand = FALSE) + scale_x_continuous(breaks = seq(-20,30,5)) + scale_y_continuous(breaks = seq(64,74,2)) + geom_rect(xmin=-3, xmax = 12, ymin = 66, ymax = 72, color = "brown1", fill = NA, size = 1) + theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))
dev.off()

x_min = 1
x_max = 70
y_min = 123
y_max = 192

x_points = seq(x_min,x_max,3)
y_points = seq(y_min,y_max,3)

x_1 = current$x[x_min,y_min]
x_2 = current$x[x_min,y_max]
x_3 = current$x[x_max,y_min]
x_4 = current$x[x_max,y_max]
y_1 = current$y[x_min,y_min]
y_2 = current$y[x_min,y_max]
y_3 = current$y[x_max,y_min]
y_4 = current$y[x_max,y_max]
d=data.frame(x=c(x_2,x_1,x_3,x_4), y=c(y_2,y_1,y_3,y_4))

pdf(file="C:/Users/henri/Documents/GitHub/Master-Thesis/Images/map_of_area_sinmod.pdf")
ggplot(data = world) + geom_sf() + coord_sf(xlim = c(0, 10), ylim = c(60, 65), expand = FALSE) + geom_polygon(data=d, mapping=aes(x=x, y=y), color = "brown1", fill = NA, size = 1) + xlab("") + ylab("") + 
  annotate(geom = "text", x = 7.8, y = 62.5, label = "Romsdalsfjorden", fontface = "italic", color = "grey22", size = 3) + annotate(geom = "text", x = 6.8, y = 60.6, label = "Hardangerfjorden", fontface = "italic", color = "grey22", size = 3) + annotate(geom = "text", x = 6, y = 61.7, label = "Nordfjord", fontface = "italic", color = "grey22", size = 3) + annotate(geom = "text", x = 6.5, y = 61.4, label = "Sognefjorden", fontface = "italic", color = "grey22", size = 3)  + theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))
dev.off()



