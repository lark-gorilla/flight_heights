## Calculation of flight height from barometric pressure for dynamic soaring seabirds

This repository accompanies the manuscript *Altimeters on albatrosses: quantifying flight heights for dynamic soaring seabirds* by Mark Miller, Sheryl Hamilton and Rohan Clarke, published in Bird Conservation International.
The code **altimeters_on_albatrosses.R** documents the manuscript analyses, while a simple worked example of the method is given below.

## Minimum worked example

The example below demonstrates the analytical approach using a simulated **5-minute burst of pressure data sampled at 1 Hz**. The simulated pressure oscillations mimic those produced by an albatross undertaking repeated dynamic soaring cycles.

### 1. Load required packages

```r
library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)
library(pracma)
library(mgcv)
```

### 2. Simulate a pressure burst

Generate 5 minutes (300 seconds) of 1 Hz pressure data with regular oscillations and random measurement noise.

```r
set.seed(123)

t <- 1:300

pres_burst <- data.frame(
  time = t,
  pressure = 101330 + 35 * sin(2 * pi * t / 15) + rnorm(300, 0, 8)
)

plot(
  pressure ~ time,
  data = pres_burst,
  type = "l",
  xlab = "Time (s)",
  ylab = "Pressure (Pa)"
)
```

### 3. Segment the burst into dynamic soaring cycles

Pressure data are first smoothed using a 3-point moving mean. Local pressure minima are then identified and used to define individual dynamic soaring cycles. A minimum peak-to-peak distance of 5 seconds is used, corresponding to the minimum dynamic soaring cycle duration reported by Schoombie et al. (2023).

```r
pres_burst$ds_seg_pressure <- NA

# Apply 3-point moving mean
pres_smth <- rollmean(
  pres_burst$pressure,
  k = 3,
  fill = NA
)

# Identify pressure minima
pres_valz <- findpeaks(
  -pres_smth,
  minpeakdistance = 5,
  nups = 2,
  ndowns = 2,
  zero = "+"
)

# Assign a unique ID to each dynamic soaring cycle
pres_burst[pres_valz[, 4] %>% sort(), ]$ds_seg_pressure <-
  seq_len(nrow(pres_valz))

pres_burst$ds_seg_pressure[1] <- 0

pres_burst <- pres_burst %>%
  fill(ds_seg_pressure, .direction = "down")
```

Plot the resulting dynamic soaring segments. Alternating black and grey lines distinguish adjacent cycles. The pressure axis is reversed so that upward movement on the plot corresponds to increasing altitude.

```r
cols <- c(
  rep(c("black", "grey"), nrow(pres_valz)),
  "black"
)

ggplot() +
  geom_line(
    aes(
      x = pres_burst$time,
      y = pres_burst$pressure,
      colour = factor(pres_burst$ds_seg_pressure)
    ),
    group = 1
  ) +
  scale_colour_manual(values = cols) +
  scale_y_reverse() +
  theme_bw() +
  theme(legend.position = "none")
```

### 4. Estimate reference pressure (`p0`)

Flight height derived from barometric pressure depends on the **reference pressure (`p0`)** assumed to represent sea level. Three scenarios are considered: **upper**, **lower**, and **central**.

#### Upper scenario

The maximum pressure recorded during the entire burst is used as `p0`. This assumes that the bird reaches (or skims) the ocean surface only once during the burst, at the point of maximum recorded pressure.

```r
pres_burst <- pres_burst %>%
  mutate(p0_mx = max(pressure))
```

#### Lower scenario

The maximum `p0` is reset separately for each dynamic soaring cycle. This assumes that the bird skims the ocean surface during **every** dynamic soaring cycle.

A temporary pressure column offset by one second is used so that a maximum occurring at the boundary between two cycles can contribute to either adjacent cycle.

```r
pres_burst$pressure1 <- pres_burst$pressure[c(1, 1:299)]

pres_burst <- pres_burst %>%
  group_by(ds_seg_pressure) %>%
  mutate(p0_ds_seg = max(pressure1, pressure)) %>%
  ungroup()

pres_burst$pressure1 <- NULL
```

#### Central scenario

A **generalised additive model (GAM)** is fitted midway between the lower and upper `p0` estimates to provide a central estimate.

```r
pres_burst$p0_diff <- pres_burst$p0_mx - pres_burst$p0_ds_seg

pres_burst$p0_gam <- fitted(
  gam(
    (p0_mx - (p0_diff / 2)) ~ s(time, k = 7),
    data = pres_burst
  )
)

# Constrain predictions to the upper and lower scenarios
pres_burst$p0_gam <- ifelse(
  pres_burst$p0_gam > pres_burst$p0_mx,
  pres_burst$p0_mx,
  pres_burst$p0_gam
)

pres_burst$p0_gam <- ifelse(
  pres_burst$p0_gam < pres_burst$p0_ds_seg,
  pres_burst$p0_ds_seg,
  pres_burst$p0_gam
)
```

### 5. Visualise the reference-pressure scenarios

The dynamic soaring segments are shown in alternating black and grey, with the three `p0` scenarios overlaid:

* **Red:** upper scenario
* **Blue:** lower scenario
* **Cyan:** central scenario

```r
cols <- c(
  rep(c("black", "grey"), nrow(pres_valz)),
  "black"
)

ggplot(pres_burst) +
  geom_line(
    aes(
      x = time,
      y = pressure,
      colour = factor(ds_seg_pressure)
    ),
    group = 1
  ) +
  scale_colour_manual(values = cols) +
  scale_y_reverse() +
  geom_line(aes(x = time, y = p0_mx), colour = 2) +
  geom_line(aes(x = time, y = p0_ds_seg), colour = 4) +
  geom_line(aes(x = time, y = p0_gam), colour = 5) +
  theme_bw() +
  theme(legend.position = "none")
```

Again, the pressure axis is reversed so that upward movement on the plot corresponds to increasing altitude.

### 6. Calculate flight height

Flight height is calculated from pressure using the barometric equation (Berberan-Santos et al. 1997):

$$
h = -\frac{kT}{mg}\ln\left(\frac{p}{p_0}\right)
$$

where:

* `h` = estimated flight height
* `p` = measured pressure
* `p0` = reference pressure
* `T` = air temperature in Kelvin
* `k` = universal gas constant
* `m` = molar mass of dry air
* `g` = gravitational acceleration

For this simulated example, air temperature is held constant at **20 °C**.

```r
k <- 8.31432
m <- 0.0289644
g <- 9.80665

example_temp <- 20

# Upper scenario
pres_burst$altitude_upper <- -(
  (k * (example_temp + 273.15)) / (m * g)
) * log(pres_burst$pressure / pres_burst$p0_mx)

# Lower scenario
pres_burst$altitude_lower <- -(
  (k * (example_temp + 273.15)) / (m * g)
) * log(pres_burst$pressure / pres_burst$p0_ds_seg)

# Central scenario
pres_burst$altitude_central <- -(
  (k * (example_temp + 273.15)) / (m * g)
) * log(pres_burst$pressure / pres_burst$p0_gam)

# Difference between upper and lower scenarios
pres_burst$altitude_diff <- -(
  (k * (example_temp + 273.15)) / (m * g)
) * log(pres_burst$p0_ds_seg / pres_burst$p0_mx)
```

### 7. View estimated flight heights

Summary statistics for flight-height estimates under the three scenarios, together with the difference between the upper and lower estimates, can be viewed using:

```r
summary(
  pres_burst[, c(
    "altitude_upper",
    "altitude_lower",
    "altitude_central",
    "altitude_diff"
  )]
)
```

---

### Interpretation

This worked example illustrates the main steps used to estimate flight height from high-frequency pressure measurements:

1. **Segment** the pressure time series into individual dynamic soaring cycles.
2. Define **upper and lower reference-pressure (`p0`) scenarios** representing alternative assumptions about when the bird approaches sea level.
3. Fit a **central `p0` estimate** between these scenarios using a GAM.
4. Convert measured pressure to **estimated flight height** using the barometric equation.
5. Use the upper and lower scenarios to quantify the **plausible range in flight-height estimates** associated with uncertainty in `p0`.
