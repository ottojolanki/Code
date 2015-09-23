##kuinka diskretisoida jatkuva muuttuja, ja laskea ryhmien keskiarvot
require(dplyr)
set.seed(123)
x <- data.frame(value = rnorm(100))
x$bins <- cut(x$value, breaks = 30)
y <- x %>% group_by(bins) %>% summarise(mean(value)) 
z <- full_join(x, y, by = "bins")
