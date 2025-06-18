# make p-values < 0.05 bold in kableExtra tables

mutate(p_value = ifelse(p.value < 0.05, 
                          ifelse(p.value < 0.001,
                                 "<b><0.001</b>", 
                                 sprintf("<b>%.3f</b>", p.value)), 
                          sprintf("%.3f", p.value)))  |> 

# escape and format are required
kable(escape = F, format = "html",