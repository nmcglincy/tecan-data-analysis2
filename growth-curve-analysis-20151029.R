library(readr)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

# TODO - make it generic so I can call it from Rscript

# skip changed to 35, how stable is it?

data = read_csv(file = "cor-titr-20151029.csv",
                skip = 35,
                n_max = 96,
                col_names = TRUE) %>%
  na.omit()
names(data)[2] = "Temp.C"
# data
# tail(data)

names(data) = make.names(names = names(data),
                         unique = TRUE,
                         allow_ = TRUE)
# data

data$Cycle.Nr. = as.numeric(data$Cycle.Nr.)

# tail(data)

data = data %>%
  select(c(1:2, 4, seq(from = 3, to = length(data), by = 2))) %>%
  gather("well", "od600", A1:H12) 

sample.info = read_csv(file = "sample-info.csv",
                       col_names = TRUE)
# sample.info

data = join(data, sample.info, by = "well")

# head(data)

cycle.data  = data %>%
  select(Cycle.Nr., Temp.C) %>%
  distinct(Cycle.Nr.)

png(file = "cycle-temp-plot.png", width = 7, height = 7, units = "in", res = 300)
plot(x = cycle.data$Cycle.Nr.,
     y = cycle.data$Temp.C,
     type = "b",
     ylim = c(29, 31),
     xlim = c(0, 96),
     xlab = "Cycle",
     ylab = "Temperature, C")
with(cycle.data, abline(lm(Temp.C ~ Cycle.Nr.), col = "blue"))
with(cycle.data, abline(line(Cycle.Nr., Temp.C), col = "red"))
dev.off()

# 
# Some statistical analysis:
# might be cool is exclude the first point, and then include some measure of variance...

cycle.temp.stats = capture.output(with(cycle.data, summary(lm(Temp.C ~ Cycle.Nr.))))
cat(cycle.temp.stats, file = "cycle-temp-stats.txt", sep = "\n", append = TRUE)
cat("Spearman's correlation coefficient", file = "cycle-temp-stats.txt", sep = "\n", append = TRUE)
cycle.temp.stats = capture.output(cor(cycle.data$Cycle.Nr. , cycle.data$Temp.C, method = "spearman"))
cat(cycle.temp.stats, file = "cycle-temp-stats.txt", sep = "\n", append = TRUE)
cat("Pearson's correlation coefficient", file = "cycle-temp-stats.txt", sep = "\n", append = TRUE)
cycle.temp.stats = capture.output(cor(cycle.data$Cycle.Nr., cycle.data$Temp.C, method = "pearson"))
cat(cycle.temp.stats, file = "cycle-temp-stats.txt", sep = "\n", append = TRUE)

# 
# head(data)
# names(data)[6] = "Row"
# names(data)[7] = "Column"

library(ggplot2)
library(ggthemes)

ggplot(data,
       aes(x = Time..ms./3600000,
           y = od600,
           colour = strain)) +
  geom_line(size = 1) +
  facet_wrap( ~ well, ncol = 12) +
  xlab("Time, hrs") +
  ylab("OD, 600 nm") +
  scale_x_continuous(breaks = c(0, 10)) +
  scale_colour_colorblind() +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=16),
        axis.text.y  = element_text(size=16),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16))
ggsave("by-well.png")

# still some noisy well, but they seem less common that before
# there's relatively few compared to the per strain per conc. replicate no, so I'll manually 
# filter them
# from IAA expt
# bad.wells = c("A2","B5","B6","C1","C2","C3","C4","C5","C6")
# iaa expt
# bad.wells = c("A1", "A5", "B5", "B6", "C1", "C6")
# iaa expt 2
# bad.wells = c("A6", "B4", "C6", "A1", "B3","C5")
bad.wells = c("A4", "B4", "E4", "E5", "F5", "F12", "D8", "D10", "F9")
good.wells = setdiff(data$well, bad.wells)

data.filt = data %>%
  filter(well %in% good.wells)

# summary(data$well)
# summary(data.filt$well)

# head(data.filt)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data.filt,
       aes(x = Time..ms./3600000,
           y = od600,
           colour = as.factor(cor.uM),
           group = well)) +
  geom_line(size = 1) +
  facet_wrap(~ strain) +
  xlab("Time, hrs") +
  ylab("OD, 600 nm") +
  scale_colour_manual(name = expression(paste("[COR], ", mu, "M", sep = "")),
                      values = cbbPalette[c(1:4,7,6)]) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        # panel.background = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=16),
        axis.text.y  = element_text(size=16),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16))
ggsave("by-strain.png")

ggplot(data.filt,
       aes(x = Time..ms./3600000,
           y = od600,
           colour = strain,
           group = well)) +
  geom_line(size = 1) +
  facet_wrap(~ cor.uM) +
  xlab("Time, hrs") +
  ylab("OD, 600 nm") +
  scale_colour_colorblind() +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=16),
        axis.text.y  = element_text(size=16),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16))
ggsave("by-conc.png")

ggplot(data.filt %>% filter(strain == "YEPD"),
       aes(x = Time..ms./3600000,
           y = od600,
           colour = as.factor(cor.uM),
           group = well)) +
  geom_line(size = 1) +
  xlab("Time, hrs") +
  ylab("OD, 600 nm") +
  scale_y_continuous(limit = c(0.1, 0.2)) +
  scale_colour_colorblind(name = expression(paste("[COR], ", mu, "M", sep = ""))) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=16),
        axis.text.y  = element_text(size=16),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16))
ggsave("yepd-wells.png")

# I don't think there is an effect of [iaa] on OD600, but I'm going to do the correction
# anyway just to have a look

data.filt.summ = data.filt %>%
  filter(strain == "YEPD") %>%
  group_by(cor.uM) %>%
  summarise(med = median(od600))

yepd.data = unlist(data.filt.summ$med)
# yepd.data

data.by.iaa.l = dlply(data.filt, .(cor.uM))
# str(data.by.iaa.l)
for (i in 1:length(data.by.iaa.l)) {
  data.by.iaa.l[[i]] = data.by.iaa.l[[i]] %>%
    mutate(c.od600 = od600 - yepd.data[i])
}
data.filt.cor = ldply(data.by.iaa.l)
# head(data.filt.cor)

# don't need the yepd data, put in strain filter here
right.strains = unique(unlist(data.filt$strain))[1:3]


data.filt.cor2 = data.filt.cor %>%
  filter(strain %in% right.strains)

ggplot(data.filt.cor2,
       aes(x = Time..ms./3600000,
           y = c.od600,
           colour = as.factor(cor.uM),
           group = well)) +
  geom_line(size = 1) +
  facet_wrap(~ strain, ncol = 2) +
  xlab("Time, hrs") +
  ylab("Background subtracted OD, 600 nm") +
  scale_colour_manual(name = expression(paste("[COR], ", mu, "M", sep = "")),
                      values = cbbPalette[c(1:4,7,6)]) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=16),
        axis.text.y  = element_text(size=16),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16))
ggsave("cor-strain.png")

ggplot(data.filt.cor2,
       aes(x = Time..ms./3600000,
           y = c.od600,
           colour = strain,
           group = well)) +
  geom_line(size = 1) +
  facet_wrap(~ cor.uM) +
  xlab("Time, hrs") +
  ylab("Background subtracted OD, 600 nm") +
  scale_colour_colorblind() +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=16),
        axis.text.y  = element_text(size=16),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16))
ggsave("cor-strain-by-conc.png")

# spline fit on corrected data

data.l2 = dlply(data.filt.cor2, .(strain, cor.uM))

growth.curve.lite2 = function(foo) {
  # calculating interesting characteristics of the growth curve, fitting a smooth.spline()
  # on an hour time-base.
  require(grofit)
  spl.fit = smooth.spline(foo$Time..ms./3600000, foo$c.od600, keep.data = FALSE)
  data.frame(fit.time = spl.fit$x, abs.fit = spl.fit$y)
} 

spl.fits.l2 = lapply(data.l2, growth.curve.lite2)

# str(spl.fits.l)

spl.fits.df2 = ldply(spl.fits.l2) %>%
  separate(.id, c("strain", "cor.uM"), sep = 7) %>%
  mutate(strain = str_sub(strain, start = 1, end = 6))

# tail(spl.fits.df2)

ggplot(data.filt.cor2, 
       aes(x = Time..ms./3600000, 
           y = c.od600,
           colour = as.factor(cor.uM))) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ strain, ncol = 2) +
  geom_line(data = spl.fits.df2, 
            aes(x = fit.time, 
                y = abs.fit, 
                fill = cor.uM), 
            size = 1.5) +
  xlab("Time, hrs") +
  ylab("Background subtracted OD, 600 nm") +
  ggtitle("smooth.spline fit by strain on pooled wells") +
  scale_colour_manual(name = expression(paste("[COR], ", mu, "M", sep = "")),
                      values = cbbPalette[c(1:4,7,6)]) +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=16),
        axis.text.y  = element_text(size=16),
        plot.title = element_text(size = 20, vjust = 1),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16))
ggsave("cor-smooth-spline-bySample.png")

