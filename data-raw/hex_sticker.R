library(tidyverse)
library(hexSticker)
library(showtext)

## Loading Google fonts (http://www.google.com/fonts)
font <- "Roboto Mono"
font_add_google(font)
## Automatically use showtext to render text for future devices
showtext_auto()

chi11_1k_expr <- read_tsv(file.path("data-raw", "chi11_1k.tsv"))

grouped <- group_by(chi11_1k_expr, sample_id, stim_type)
one_group <- group_split(grouped)[[1]]

state_marker <- "pSTAT3"
one_marker <- one_group[state_marker]

p <- ggplot(one_marker, aes(x = !!ensym(state_marker))) +
  geom_density(color="#0F3B13", fill="#486144")
p <- p + theme_void() + theme_transparent()

s <- sticker(p, package="HDStIM", p_size=20.5, p_family = font, p_x = 1, p_y = 1.45,
             s_x=1, s_y=.8, s_width=1, s_height=0.9,
             h_fill="#E47833", h_color="#f39c12", h_size = 1, dpi = 300,
             filename=file.path("man","figures","sticker.png"))

# https://www.december.com/html/spec/color1.html

plot(s)

