# Samuel Moijueh
# February 2019

# load required packages
library(pacman)
pacman::p_load(survival, survminer, dplyr, ggplot2, gridExtra, DT, corrplot, broom)

# import PBC dataset, encode missing values as 'NA'
d <- read.table(file = "https://www4.stat.ncsu.edu/~boos/var.select/pbc.dat.txt",
                       na.strings = c('.', 'NA'))

# set column names
colnames(d) <- c('id','futime','status','drug','age','sex','ascites',
                        'hepato','spiders','edema','bili','chol','albumin',
                        'copper','alk_phos','sgot','trig','platelet',
                        'protime','stage')

# convert the time in days to years for simplicity
d$age = floor(d$age / 365)

# remove all rows with missing values from the dataframe
d <- d[complete.cases(d), ]

# recode status
d$status = ifelse(d$status == 2, 1, 0)
table(d$status)/length(d$status)

# Create a histogram
ggplot(d, aes(x=age)) +
    geom_histogram(binwidth=2, colour="black", fill="pink") + ggtitle("Distribution of Age") + xlab("Age (years)") +  ylab("Frequency") +   ggtitle("Histogram of Age") +
  theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks = seq(0, 80, by = 3))


# correlation matrix

cor.mtest <- function(mat) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j])
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of correlation statistic
corr<-cor(d)

# matrix of the p-value of the correlation
p.mat <- cor.mtest(d)

#corrplot(corr, method = "circle")

# Specialized the insignificant value according to the significant level
corrplot(corr, type="upper", order="hclust",
         p.mat = p.mat, sig.level = 0.01)  # 0.01 is more stringent than 0.05.
# small p-value reject null hypothesis.

# density distribution
xlabels <- list("bili" = "Serum Bilirubin (mg/dl)",
                "chol" = "Serum Cholesterol in (mg/dl)",
                "albumin" = "Albumin (gm/dl)",
                "copper" = "Urine Copper (ug/day))",
                "alk_phos" = "Alkaline Phosphatase (U/liter))",
                "sgot" = "SGOT (U/ml))",
                "trig" = "triglicerides in mg/dl platelet = platelets per cubic ml/1000",
                "protime" = "prothrombin time in seconds")

density_plot<-function(x, log = FALSE){
  if(log){
    return(ggplot(d, aes(x=log(d[[x]]))) + geom_density() + xlab(paste0("log(d$",x,")")) +  ylab("Density") + ggtitle(paste0("log(d$", x,")")) + theme(plot.title = element_text(hjust = 0.5)))
  }
  return(ggplot(d, aes(x=d[[x]])) + geom_density() + xlab(xlabels[[x]]) +  ylab("Density") + ggtitle(paste0("Density Plot of d$", x)) + theme(plot.title = element_text(hjust = 0.5)))
}

l <- ls(xlabels)

# Arrange and display the plots into a 2x4 grid. xlab(paste0("log( d$", x,")"))
grid.arrange(density_plot(l[1]),
             density_plot(l[2]),
             density_plot(l[3]),
             density_plot(l[4]), ncol=2)

grid.arrange(density_plot(l[5]),
             density_plot(l[6]),
             density_plot(l[7]),
             density_plot(l[8]), ncol=2)

# apply log transformation on left skewed data
# Arrange and display the plots into a 2x1 grid
grid.arrange(density_plot(l[2], 1),
             density_plot(l[3], 1),
             density_plot(l[4], 1),
             density_plot(l[5], 1),
             density_plot(l[6], 1),
             density_plot(l[7], 1),
             density_plot(l[8], 1), ncol=4)

# apply log transformation -- 7 variables
d$alk_phos <- log(d$alk_phos)
d$bili <- log(d$bili)
d$chol <- log(d$chol)
d$copper <- log(d$copper)
d$protime <- log(d$protime)
d$sgot <- log(d$sgot)
d$trig <- log(d$trig)

# change data labels -- 7 variables
d$drug <- factor(d$drug,
                     levels = c("1", "2"),
                     labels = c("D-penicillamine", "placebo"))

d$sex <- factor(d$sex,
                     levels = c("0", "1"),
                     labels = c("male", "female"))

d$ascites <- factor(d$ascites,
                    levels = c("0", "1"),
                    labels = c("no", "yes"))

d$hepato <- factor(d$hepato,
                    levels = c("0", "1"),
                    labels = c("no", "yes"))

d$spiders <- factor(d$spider,
                    levels = c("0", "1"),
                    labels = c("no", "yes"))

d$edema <- factor(d$edema,
                    levels = c("0", "0.5", "1"),
                    labels = c("no edema and no diuretic therapy for edema",
                               "edema present without diuretics, or edema resolved by diuretics",
                               "edema despite diuretic therapy"))

d$stage <- factor(d$stage,
                      levels = c("1", "2", "3", "4"),
                      labels = c("one", "two", "three", "four"))

## Estimating Survival Time -- Life Table
# create the survival object. the event is denoted by 1
fusurv <- Surv(d$futime, d$status)
fusurv

# Fitting the Kaplan-Meier Estimate
fit_KM = survfit(fusurv ~ drug, type="kaplan-meier", data = d)

# Life Table
summary(fit_KM)

# Kaplan-Meier Survival Curve and Log-Rank Hypothesis Test
ggsurv<-ggsurvplot(fit_KM, data=d,
                   pval = TRUE,
                   conf.int = TRUE,
                   ggtheme = theme_bw(),
                   surv.median.line = "v",
                   legend.labs = c("D-penicillamine", "Placebo"),
                   title = "Survival curves",
                   subtitle = "Based on Kaplan-Meier estimates",
                   xlab = "Time in days")
                   #palette = c("#E7B800", "#2E9FDF")) # custom color palettes.)

# Labels for Survival Curves (plot)
ggsurv$plot <- ggsurv$plot + labs(
  title    = "\tSurvival curves",
  subtitle = "\tBased on Kaplan-Meier estimates",
  caption = NULL
  ) + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.subtitle = element_text(hjust = 0.5))

ggsurv <- ggpar(
  ggsurv,
  font.title    = c(16,"darkblue"),
  font.subtitle = c(15, "purple")
)

ggsurv

# Estimate x-day survival rate
# 2256 days is about 6 years and four months
summary(fit_KM, times = 2256)

# Cox Model Assumptions
res.cox <- coxph(fusurv ~ drug + age + sex + ascites + hepato +
                   spiders + edema + bili + chol + albumin + copper + alk_phos +
                   sgot + trig + protime + stage , data = d)

test.ph <- cox.zph(res.cox)
test.ph

ggcoxfunctional(fusurv ~ age + bili + chol + albumin + copper +
                  alk_phos + sgot + trig + protime, data = d)

# Cox Proportional-Hazards Model

# Fit a Cox proportional-hazards model
res.cox <- coxph(fusurv ~ drug + age + sex + ascites + hepato +
                   spiders + edema + bili + chol + albumin + copper + alk_phos +
                   sgot + trig + protime + stage , data = d)

broom::tidy(res.cox)

# Visualizing the Hazard Ratios via Forest Plot

ggforest(res.cox, data = d)
