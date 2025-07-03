
test <- d04_sel1 |> 
  # mutate(poradie_vysetrenia = str_remove(poradie_vysetrenia, "M") |> as.numeric()) |>
  dplyr::select(-odpoved_na_terapii_m0_vs_m6) |> 
  pivot_longer(cols = any_of(var_dep_01),
               names_to = "var_dep_name",
               values_to = "var_dep_value") |> 
  pivot_longer(cols = any_of(var_indep_01),
               names_to = "var_indep_name",
               values_to = "var_indep_value") |> 
  group_by(var_dep_name, var_indep_name) |> 
  nest() |> 
  mutate(mod_01 = map(data, ~possLMER(data = .x,
                                      family = lognormal(link = "log"),
                                      var_dep_value ~ poradie_vysetrenia + var_indep_value + (1|projekt_id))),
         tidier_01 = map(mod_01,
                         broom.mixed::tidy,
                         conf.int   = TRUE,
                         conf.level = 0.95),
         mod_02 = map(data, ~possLMER(data = .x,
                                      family = lognormal(link = "log"),
                                      var_dep_value ~ poradie_vysetrenia + var_indep_value + (poradie_vysetrenia|projekt_id))),
         tidier_02 = map(mod_02,
                         broom.mixed::tidy,
                         conf.int   = TRUE,
                         conf.level = 0.95)
         )




test$mod_01[[1]] |> check_model()
test$mod_02[[1]] |> check_model()

test$mod_01[[1]] |> tab_model()
test$mod_02[[1]] |> tab_model()


data_test <- test$data[[1]]



#####################



# 1) Nainstalujte a načtěte brms
# install.packages("brms")
library(brms)

# 2) Vytvořte indikátor cenzury
#    hodnoty ≥ 80 budou označeny jako "right" (pravá cenzura), ostatní "none"
data_test <- data_test |> 
  mutate(
    cens = ifelse(var_dep_value >= 80, "right", "none"),
    # doporučuji také škálovat kovariáty
    poradie_sc = as.numeric(scale(poradie_vysetrenia)),
    indep_sc   = as.numeric(scale(var_indep_value))
  )

# 3) Specifikujte model pomocí bf(), kde | cens(cens) říká, co je cenzurované
tobit_bf <- bf(
  var_dep_value | cens(cens) ~ 
    poradie_sc + indep_sc + 
    (poradie_sc | projekt_id)
)

# 4) (Volitelné) Definujte priory
priors <- c(
  prior(normal(0, 1), class = "b"),         # slope priors
  prior(exponential(1), class = "sigma"),   # prior na σ
  prior(lkj(2), class = "cor")              # prior na korelaci random efektů
)

detach("package:conflicted", unload = TRUE) # I recommend detaching otherwise you spend a lot of time solving the conflicts
# 5) Fitting modelu
fit_brms <- brm(
  formula = tobit_bf,
  family  = gaussian(),       # cens(cens) v bf() aktivuje Tobit likelihood
  data    = data_test,
  prior   = priors,
  chains  = 4,
  cores   = ncores,
  iter    = 4000,
  control = list(adapt_delta = 0.95),
  sample = TRUE,        # implicitně ano, nemusíte psát
  algorithm = "sampling" # implicitně MCMC
)

# 6) Shrnutí výsledků
print(summary(fit_brms))
plot(fit_brms)
pp_check(fit_brms, nsamples = 100)





library(censReg)
library(plm)

# 1) Převod na panel data frame (list) s id a „časem“
pData <- pdata.frame(
  data_test,
  index = c("projekt_id", "poradie_vysetrenia")
)

# 2) poradie_vysetrenia se v rámci pdata.frame stane faktorem,
#    takže ho musíme zpátky na numerickou:
pData$poradie_num <- as.numeric(as.character(pData$poradie_vysetrenia))

# 3) Fit censReg – pravá cenzura na 80
fit_cr <- censReg(
  formula = var_dep_value ~ poradie_vysetrenia + var_indep_value,
  data    = pData,
  left    = 0,
  right   = 80,
  method  = "BHHH"
)

# 4) Výsledky
summary(fit_cr)





library(brms)
library(bayestestR)

data_test <- data_test |> 
  mutate(
    cens = ifelse(var_dep_value >= 80, "right", "none"),
    # doporučuji také škálovat kovariáty
    poradie_sc = as.numeric(scale(poradie_vysetrenia)),
    indep_sc   = as.numeric(scale(var_indep_value))
  )


data_test <- test$data[[1]]
data_test <- data_test %>%
  mutate(
    y1        = var_dep_value,
    y2        = ifelse(var_dep_value >= 80, 80, var_dep_value),
    cens_type = factor(
      ifelse(var_dep_value >= 80, "right", "none"),
      levels = c("none","left","right","interval")
    )
  )

mod1 <- brm(y1 | cens(cens_type, y2) ~ poradie_vysetrenia + var_indep_value + (1|projekt_id),
            data = data_test,
            family = gaussian(), chains=4, cores=ncores)

summary(mod1)
plot(mod1)
conditional_effects(mod1, "var_indep_value")

plot(conditional_effects(mod1, effects = "var_indep_value",
                         conditions = conditions, re_formula = NULL),
     points = TRUE, rug = TRUE)



data_test |> 
  ggplot(aes(var_indep_value, exp(var_dep_value))) +
  geom_point() +
  stat_smooth(method = "loess")



glmmTMB(var_indep_value ~ poradie_vysetrenia + var_indep_value + (1|projekt_id), 
        family = tobit(Upper = 80), data = data_test)

mixed_model(family = censored.normal(censored.normal(right = 80)),
            fixed = var_indep_value ~ poradie_vysetrenia + var_indep_value, 
            random = ~1|projekt_id,
            data = data_test
            )





