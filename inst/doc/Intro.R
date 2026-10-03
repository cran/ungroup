## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.align = "center",
  dpi       = 96,
  out.width = "90%"
)

## ----load---------------------------------------------------------------------
library(ungroup)

## ----overview-diagram, echo = FALSE, fig.height = 3.2, fig.width = 8----------
op <- par(mar = c(0.2, 0.2, 0.2, 0.2))
plot.new()
# Room on the left for the labels, or they are clipped at the margin.
plot.window(xlim = c(-3.4, 10.4), ylim = c(0, 3.6))

# Observed bins: unequal widths, the last one wide and open-ended.
bin_width <- c(1, 1, 3, 5, 10, 10) / 3   # scaled so the total span is 10
bin_counts <- c(294, 66, 32, 170, 284, 998)
left <- cumsum(c(0, head(bin_width, -1)))

for (i in seq_along(bin_counts)) {
  rect(left[i], 2.4, left[i] + bin_width[i], 3.0,
       col = "grey80", border = "white")
  text(left[i] + bin_width[i] / 2, 3.15, bin_counts[i], cex = 0.7)
}
text(-0.25, 2.7, "observed counts", adj = 1, cex = 0.8)

text(5, 2.05, expression(y %~% Poisson(C * gamma)), cex = 1.1)
arrows(5, 1.9, 5, 1.55, length = 0.1)

# The fine grid the model estimates, spanning exactly the same range.
n_cell <- 40
cw <- 10 / n_cell
for (i in seq_len(n_cell)) {
  rect((i - 1) * cw, 0.9, i * cw, 1.5,
       col = if (i %% 2) "steelblue" else "lightsteelblue", border = "white")
}
text(-0.25, 1.2, "estimated", adj = 1, cex = 0.8)

text(5, 0.45, "fine grid, one value per cell", cex = 0.8)
text(5, 0.1, "C aggregates the fine grid back onto the coarse bins",
     cex = 0.75, col = "grey30")
par(op)

## ----data---------------------------------------------------------------------
# x: start of each input interval. The last interval runs [85, 85 + nlast).
x <- c(0, 1, seq(5, 85, by = 5))

# y: deaths in each interval
y <- c(294, 66, 32, 44, 170, 284, 287, 293, 361, 600, 998,
       1572, 2529, 4637, 6161, 7369, 10481, 15293, 39016)

# nlast: width of the open final interval, so [85, 111)
nlast <- 26

## ----quickstart---------------------------------------------------------------
M1 <- pclm(x = x, y = y, nlast = nlast)
M1

## ----quickstart-plot----------------------------------------------------------
plot(M1, xlab = "Age, x", ylab = "Deaths")

## ----object-------------------------------------------------------------------
names(M1)

## ----object-summary-----------------------------------------------------------
summary(M1)

## ----residuals----------------------------------------------------------------
head(residuals(M1), 5)
max(abs(residuals(M1)))

## ----omega--------------------------------------------------------------------
M_omega <- pclm(x = x, y = y, omega = 111, control = list(lambda = 100))
M_nlast <- pclm(x = x, y = y, nlast = 26,  control = list(lambda = 100))
all.equal(unname(fitted(M_omega)), unname(fitted(M_nlast)))

## ----omega-errors, error = TRUE-----------------------------------------------
try({
pclm(x = x, y = y, nlast = 26, omega = 111)  # both given
pclm(x = x, y = y)                           # neither given
pclm(x = x, y = y, omega = 85)               # not past max(x)
})

## ----outstep------------------------------------------------------------------
M2 <- pclm(x = x, y = y, nlast = nlast, out.step = 0.5)
length(fitted(M2))
head(names(fitted(M2)), 4)

## ----outstep-mass-------------------------------------------------------------
c(same_mass = all.equal(sum(fitted(M2)), sum(y)))

## ----outstep-warning----------------------------------------------------------
M2b <- pclm(x = x, y = y, nlast = nlast, out.step = 0.32,
            control = list(lambda = 100))

## ----outstep-suggest----------------------------------------------------------
suggest.valid.out.step(max(x) + nlast - min(x))

## ----e0-setup-----------------------------------------------------------------
# Average years lived, from a vector of deaths by single year of age.
e0 <- function(dx) {
  n <- length(dx)
  l <- rev(cumsum(rev(dx)))
  l <- l / l[1]
  L <- c((l[-1] + l[-n]) / 2, l[n])
  sum(L) / l[1]
}

# Aggregate single-age deaths into the same coarse bins used above.
grp <- rep(x, c(diff(x), nlast))
bin_deaths <- function(j) {
  as.numeric(tapply(as.numeric(ungroup.data$Dx[, j]), grp, sum))
}

## ----e0-1980------------------------------------------------------------------
truth <- as.numeric(ungroup.data$Dx[, 1])
coarse <- bin_deaths(1)
widths <- c(diff(x), nlast)

uniform <- rep(coarse / widths, times = widths)
fit1980 <- pclm(x, coarse, nlast, control = list(lambda = 100))

round(c(truth     = e0(truth),
        uniform   = e0(uniform),
        pclm      = e0(unname(fitted(fit1980)))), 3)

## ----e0-all-------------------------------------------------------------------
errors <- t(vapply(1:35, function(j) {
  truth  <- as.numeric(ungroup.data$Dx[, j])
  coarse <- bin_deaths(j)
  uniform <- rep(coarse / widths, times = widths)
  fit <- pclm(x, coarse, nlast, control = list(lambda = 100))
  c(uniform = e0(uniform) - e0(truth),
    pclm    = e0(unname(fitted(fit))) - e0(truth))
}, numeric(2)))

round(apply(abs(errors), 2, function(z) c(mean = mean(z), max = max(z))), 3)

## ----e0-plot------------------------------------------------------------------
plot(1980:2014, errors[, "uniform"], type = "b", pch = 19, col = "grey60",
     ylim = range(errors) * 1.1,
     xlab = "Year", ylab = "Error in e0, years")
lines(1980:2014, errors[, "pclm"], type = "b", pch = 19, col = 2)
abline(h = 0, lty = 3)
legend("topleft", legend = c("uniform spread", "pclm"),
       col = c("grey60", 2), lty = 1, pch = 19, bty = "n")

## ----faithful-----------------------------------------------------------------
faithful_counts <- hist(datasets::faithful$eruptions,
                        breaks = seq(1.5, 5.5, by = 1),
                        plot   = FALSE)$counts
faithful_counts

Mf <- pclm(x = 1.5:4.5, y = faithful_counts, nlast = 1,
           out.step = 0.1)
Mf$smoothPar[1]

## ----faithful-plot------------------------------------------------------------
plot(Mf, xlab = "Eruption length, minutes", ylab = "Eruptions")

## ----faithful-lambda----------------------------------------------------------
valley <- function(L) {
  fv <- as.numeric(fitted(pclm(1.5:4.5, faithful_counts, 1, out.step = 0.1,
                               control = list(lambda = L))))
  round(c(valley = min(fv[11:20]), peak = max(fv)), 2)
}
rbind(lambda_1   = valley(1),
      lambda_auto = valley(Mf$smoothPar[1]),
      lambda_1e6  = valley(1e6))

## ----faithful-bic-------------------------------------------------------------
sapply(c(1, Mf$smoothPar[1], 1e6), function(L) {
  M <- pclm(1.5:4.5, faithful_counts, 1, out.step = 0.1,
            control = list(lambda = L))
  round(c(BIC = BIC(M), AIC = AIC(M)), 2)
})

## ----lambda-bound-------------------------------------------------------------
M_auto <- pclm(x, y, nlast, control = list(lambda = NA))
M_wide <- pclm(x, y, nlast,
               control = list(lambda = NA, int.lambda = c(1e-4, 1e5)))
c(default_search = M_auto$smoothPar[1], wider_search = M_wide$smoothPar[1])

## ----lambda-bic---------------------------------------------------------------
M_bic <- pclm(x, y, nlast, control = list(lambda = NA, opt.method = "BIC"))
M_aic <- pclm(x, y, nlast, control = list(lambda = NA, opt.method = "AIC"))
c(bic_choice = M_bic$smoothPar[1],
  aic_choice = M_aic$smoothPar[1])

## ----offset-------------------------------------------------------------------
Ex <- c(114, 440, 509, 492, 628, 618, 576, 580, 634, 657,
        631, 584, 573, 619, 530, 384, 303, 245, 249) * 1000

M3 <- pclm(x = x, y = y, nlast = nlast, offset = Ex)
fitted(M3)[1:5]

## ----offset-plot--------------------------------------------------------------
plot(M3, type = "s", xlab = "Age, x", ylab = "m(x), log scale")

## ----offset-two-forms---------------------------------------------------------
Ex_fine <- fitted(pclm(x = x, y = Ex, nlast = nlast))   # 111 values

M_coarse <- pclm(x = x, y = y, nlast = nlast, offset = Ex)
M_fine   <- pclm(x = x, y = y, nlast = nlast, offset = Ex_fine)

# Same rates to within a few percent through the bulk of the distribution,
# and further apart in the extreme tail where the counts are thin.
round(range(fitted(M_fine) / fitted(M_coarse)), 3)

## ----zeros--------------------------------------------------------------------
small <- c(0, 0, 1, 0, 2, 1, 0, 3, 2, 1)
M_small <- pclm(x = 0:9, y = small, nlast = 1,
                control = list(lambda = 10))

## ----zeros-scaled-------------------------------------------------------------
M_scaled <- pclm(x = 0:9, y = small * 100, nlast = 1,
                 control = list(lambda = 10))
c(total_in   = sum(small * 100),
  total_out  = sum(fitted(M_scaled)),
  head_fit   = round(head(fitted(M_scaled), 3), 1))

## ----ci-names-----------------------------------------------------------------
names(M1$ci)

## ----ci-crossing--------------------------------------------------------------
diff <- fitted(M1) - M1$ci$lower
c(totals_equal = isTRUE(all.equal(sum(M1$ci$lower), sum(y))),
  first_age_above = names(fitted(M1))[min(which(diff < 0))])

## ----ci-table-----------------------------------------------------------------
i <- c(1, 30, 70, 100, 111)
data.frame(
  bin       = names(fitted(M1))[i],
  fitted    = round(fitted(M1)[i], 1),
  conf_lo   = round(M1$ci$conf_lower[i], 1),
  conf_up   = round(M1$ci$conf_upper[i], 1),
  scen_lo   = round(M1$ci$lower[i], 1),
  scen_up   = round(M1$ci$upper[i], 1)
)

## ----ci-plot------------------------------------------------------------------
lo <- M1$ci$conf_lower
up <- M1$ci$conf_upper
f  <- fitted(M1)
age <- seq(0, 110, length.out = length(f))

plot(age, f, type = "l", lwd = 2, ylim = c(0, max(up) * 1.05),
     xlab = "Age, x", ylab = "Deaths")
polygon(c(age, rev(age)), c(lo, rev(up)),
        col = adjustcolor("steelblue", 0.25), border = NA)
lines(age, f, lwd = 2)
lines(age, M1$ci$lower, lwd = 2, lty = 2, col = 2)
lines(age, M1$ci$upper, lwd = 2, lty = 2, col = 2)
legend("topright", bty = "n", lty = c(1, 1, 2), lwd = 2,
       col = c(1, "steelblue", 2),
       legend = c("fitted", "pointwise 95%", "mass-preserving scenarios"))

## ----ragged-------------------------------------------------------------------
grp2 <- rep(x, c(diff(x), nlast))
years <- 1:12
y2d <- aggregate(ungroup.data$Dx[, years], by = list(grp2), FUN = "sum")[, -1]

# The top of the surface is unobserved in the early years: the oldest ages
# were folded into the open interval at a lower ceiling then.
y_ragged <- y2d
y_ragged[17:19, 1:4] <- NA
y_ragged[c(1:3, 16:19), 1:5]

## ----ragged-fail, error = TRUE------------------------------------------------
try({
pclm2D(x = x, y = y_ragged, nlast = nlast, verbose = FALSE)
})

## ----ragged-omit--------------------------------------------------------------
P_ragged <- pclm2D(x = x, y = y_ragged, nlast = nlast,
                   na.action = "omit", verbose = FALSE,
                   control = list(lambda = c(1, 1), kr = 5))
dim(fitted(P_ragged))
all(is.finite(fitted(P_ragged)))

## ----missing-year-------------------------------------------------------------
Ex2d <- aggregate(ungroup.data$Ex[, years], by = list(grp2), FUN = "sum")[, -1]
Ex2d[, c(2, 3, 5, 6, 8, 9, 11, 12)] <- NA

## ----missing-year-fit---------------------------------------------------------
P_missing <- pclm2D(x = x, y = y2d, nlast = nlast, offset = Ex2d,
                    na.action = "omit", verbose = FALSE,
                    control = list(lambda = c(1, 1), kr = 5))
dim(fitted(P_missing))
all(is.finite(fitted(P_missing)))

## ----omit-mass----------------------------------------------------------------
c(
  observed_cells   = sum(y_ragged, na.rm = TRUE),
  fitted_all_cells = sum(fitted(P_ragged))
)

## ----two-dimensional----------------------------------------------------------
years10 <- 1:10
y10  <- aggregate(ungroup.data$Dx[, years10], by = list(grp2), FUN = "sum")[, -1]
Ex10 <- aggregate(ungroup.data$Ex[, years10], by = list(grp2), FUN = "sum")[, -1]
dim(y10)

P_counts <- pclm2D(x = x, y = y10, nlast = nlast, verbose = FALSE,
                   control = list(lambda = c(1, 1), kr = 5))
dim(fitted(P_counts))

P_rates <- pclm2D(x = x, y = y10, nlast = nlast, offset = Ex10,
                  verbose = FALSE, control = list(lambda = c(1, 1), kr = 5))
summary(P_rates)

## ----kr-error, error = TRUE---------------------------------------------------
try({
pclm2D(x = x, y = y10[, 1:5], nlast = nlast, verbose = FALSE,
       control = list(lambda = c(1, 1)))
})

## ----kr-cost------------------------------------------------------------------
basis_size <- vapply(c(2, 3, 5, 7), function(k) {
  P <- pclm2D(x = x, y = y10, nlast = nlast, verbose = FALSE,
              control = list(lambda = c(1, 1), kr = k))
  ncol(P$deep$B)
}, numeric(1))
setNames(basis_size, paste0("kr=", c(2, 3, 5, 7)))

## ----two-dimensional-plots, fig.height = 4.8----------------------------------
plot(P_counts, xlab = "Age", ylab = "Year", zlab = "Deaths")

## ----two-dimensional-rates-plot, fig.height = 4.8-----------------------------
plot(P_rates, xlab = "Age", ylab = "Year", zlab = "log m(x)")

## ----two-dimensional-observed, fig.height = 4.8-------------------------------
plot(P_counts, type = "observed", xlab = "Age", ylab = "Year",
     zlab = "Deaths per year of age")

## ----two-dimensional-heatmap, fig.height = 4.4--------------------------------
Z <- fitted(P_rates)
image(x = as.numeric(sub("^\\[([0-9.]+),.*$", "\\1", rownames(Z))),
      y = years10, z = log(Z), col = hcl.colors(64, "YlOrBr", rev = TRUE),
      xlab = "Age, x", ylab = "Year")
contour(x = as.numeric(sub("^\\[([0-9.]+),.*$", "\\1", rownames(Z))),
        y = years10, z = log(Z), add = TRUE, col = "grey30", labcex = 0.7)

## ----session------------------------------------------------------------------
sessionInfo()

