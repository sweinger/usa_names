library(vegan)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)

# read names data from ssa's website: https://www.ssa.gov/oact/babynames/names.zip and import into data frame
temp <- tempfile()
download.file("https://www.ssa.gov/oact/babynames/names.zip", temp)
extractdir <- paste0(tempdir(),"/names")
unzip(temp, exdir = extractdir)
ssa_files <- list.files(path = extractdir, pattern = ".txt", full.names = TRUE)

baby_names <- do.call(rbind, lapply(ssa_files, function(x) {
  year <- as.integer(regmatches(x, regexec("yob(.*).txt", x))[[1]][2])
  cbind(year, read.csv(x, header = FALSE))
}))
names(baby_names) <- c("year","name","gender","count")
baby_names <- baby_names[]

unlink(temp)
unlink(extractdir, recursive = TRUE)

# step 1: compute raw diversity metrics for each year and gender and save into data frame for plotting
diversity.metrics <- data.frame(year=c(), type=c(), esn=c())

years <- unique(baby_names$year)
years <- years[order(years)]

for (y in years) {
  diversity.metrics <- rbind(diversity.metrics, data.frame(year=y,type="male",esn=exp(diversity(baby_names[baby_names$year==y&baby_names$gender=="M",]$count))))
  diversity.metrics <- rbind(diversity.metrics, data.frame(year=y,type="female",esn=exp(diversity(baby_names[baby_names$year==y&baby_names$gender=="F",]$count))))
}

p <- ggplot(subset(diversity.metrics, year >= 1937)) + 
  geom_line(aes(x = year, y = esn, group = type), size=0.8) + 
  geom_text(data = subset(diversity.metrics, year == "2016"), aes(label = type, x = Inf, y = esn, fontface = "bold"), hjust = -.1) +
  ylab(label="Effective number of species") + 
  theme(plot.margin = unit(c(1,3,1,1), "lines"), axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(face = "bold"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.title.y = element_text(angle=90,vjust =2)) 

# Code to turn off clipping
gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.newpage()
grid.draw(gt)

good_turing_counts <- function(x) {
  n <- tabulate(x + 1L)
  n0 <- n[1]
  n <- n[-1]
  max.x <- length(n)
  r <- seq.int(from = 1L, to = max.x)
  n.pos <- n[n > 0]
  r.pos <- r[n > 0]
  l <- length(n.pos)
  q <- diff(c(0L, r.pos, 2L * r.pos[l] - r.pos[l - 1]), lag = 2)/2
  z <- n.pos/q
  t_data <- data.frame(z=z,r=r.pos)
  fit <- lm(log(z) ~ log(r), data = t_data)
  return(exp(predict(fit, newdata = data.frame(r=1:4))))  
}

good_turing_plot <- function(x) {
  n <- tabulate(x + 1L)
  n0 <- n[1]
  n <- n[-1]
  max.x <- length(n)
  r <- seq.int(from = 1L, to = max.x)
  n.pos <- n[n > 0]
  r.pos <- r[n > 0]
  l <- length(n.pos)
  q <- diff(c(0L, r.pos, 2L * r.pos[l] - r.pos[l - 1]), lag = 2)/2
  z <- n.pos/q
  t_data <- data.frame(z=z,r=r.pos)
  m <- lm(log(z) ~ log(r), data = t_data)
  ggplot(t_data, aes(x = log(r), y= log(z))) + 
    geom_point() +
    geom_point(data = data.frame(z=predict(m, newdata = data.frame(r=(1:4))), r=log(1:4)),
               aes(x = r, y = z), colour = "blue") +
    stat_smooth(method = "lm", col = "red", se = FALSE) +
    theme(plot.margin = unit(c(1,3,1,1), "lines"),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour="#f0f0f0"),
          axis.text = element_text(face = "bold"),
          axis.line = element_line(colour="black"),
          axis.ticks = element_line(),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.y = element_text(angle=90,vjust =2))     
}

# log log plot for 2016 male names
good_turing_plot(baby_names[baby_names$year == 2016 & baby_names$gender == "M",]$count)
good_turing_counts(baby_names[baby_names$year == 2016 & baby_names$gender == "M",]$count)

# uses the good turing frequency estimation procedure to estimate the counts for names with less than 5 observations in a year
estimated_esn <- function(year, gender) {
  name_data <- baby_names[baby_names$year==year&baby_names$gender==gender,]
  estimated_counts <- good_turing_counts(name_data$count)
  estimated_total <- sum(name_data$count) + sum(estimated_counts*(1:4))
  proportions <- name_data$count/estimated_total
  estimated_proportions <- (1:4)/estimated_total
  esn <- exp(-sum(proportions * log(proportions)) - sum(estimated_counts*estimated_proportions*log(estimated_proportions)))
  return(list(esn=esn, estimated_total=estimated_total))
}

diversity.metrics.estimated <- data.frame(year=c(), type=c(), esn=c(), modeled_total=c())

for (y in years) {
  male_estimates <- estimated_esn(year = y, gender = "M")
  female_estimates <- estimated_esn(year = y, gender = "F")
  diversity.metrics.estimated <- rbind(diversity.metrics.estimated, data.frame(year=y,type="male",esn=male_estimates$esn, modeled_total=male_estimates$estimated_total))
  diversity.metrics.estimated <- rbind(diversity.metrics.estimated, data.frame(year=y,type="female",esn=female_estimates$esn, modeled_total=female_estimates$estimated_total))
}

# plot that compares the two methodologies: raw vs. adjusted counts
p <- ggplot(subset(diversity.metrics, year >= 1937)) + 
  geom_line(aes(x = year, y = esn, group = type), size=0.8) + 
  geom_text(data = subset(diversity.metrics, year == "2016"), aes(label = type, x = Inf, y = esn, fontface = "bold"), hjust = -.1) +
  ylab(label="Effective number of species") + 
  scale_y_continuous(limits = c(0,3000)) +
  labs(title = "Raw data") +
  theme(plot.title = element_text(face = "bold",
                                  size = rel(1.2), hjust = 0.5),
        plot.margin = unit(c(1,3,1,1), "lines"), axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(face = "bold"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.title.y = element_text(angle=90,vjust =2))
    

q <- ggplot(subset(diversity.metrics.estimated, year >= 1937)) + 
  geom_line(aes(x = year, y = esn, group = type), size=0.8) + 
  geom_text(data = subset(diversity.metrics.estimated, year == "2016"), aes(label = type, x = Inf, y = esn, fontface = "bold"), hjust = -.1) +
  scale_y_continuous(limits = c(0,3000)) +
  labs(title = "Adjusted data") +
  theme(plot.title = element_text(face = "bold",
                                  size = rel(1.2), hjust = 0.5),
        plot.margin = unit(c(1,3,1,1), "lines"), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(face = "bold"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        axis.title = element_text(face = "bold",size = rel(1)))
    #theme(plot.margin = unit(c(1,3,1,1), "lines"), axis.title.x=element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +

# Code to turn off clipping
gtp <- ggplotGrob(p)
gtp$layout$clip[gtp$layout$name == "panel"] <- "off"
gtq <- ggplotGrob(q)
gtq$layout$clip[gtq$layout$name == "panel"] <- "off"

# Code to turn off clipping
grid.newpage()
grid.arrange(gtp,gtq, nrow = 1)

# find ratio of diversity of girl names vs. boy names
gender.diversity <- reshape(diversity.metrics.estimated, idvar = "year", timevar = "type", direction = "wide")
gender.diversity$ratio <- gender.diversity$esn.female/gender.diversity$esn.male

p <- ggplot(subset(gender.diversity, year >= 1937)) + 
  geom_line(aes(x = year, y = ratio), size=0.8) + 
  ylab(label="Ratio of boy to girl names diversity") + 
  theme(plot.title = element_text(face = "bold",
                                  size = rel(1.2), hjust = 0.5),
        plot.margin = unit(c(1,3,1,1), "lines"), axis.title.x=element_blank(),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(face = "bold"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.title.y = element_text(angle=90,vjust =2))
p
