# Working models ----
library(censReg)
library(plm)

# 1) Převod na panel data frame (list) s id a „časem“
pData <- pdata.frame(
  data_test,
  index = c("projekt_id", "poradie_vysetrenia")
)

# 2) poradie_vysetrenia se v rámci pdata.frame stane faktorem,
#    takže ho musíme zpátky na numerickou:
# pData$poradie_num <- as.numeric(as.character(pData$poradie_vysetrenia))

# 3) Fit censReg – pravá cenzura na 80
fit_cr <- censReg(
  formula = var_dep_value ~ poradie_vysetrenia + log(var_indep_value),
  data    = pData,
  left    = 0,
  right   = 80,
  method  = "BHHH"
)

# 4) Výsledky
summary(fit_cr)

plot(fit_cr)


fit_cr$