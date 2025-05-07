
tidy_RD <- function(data, f1, f2, f3, disease_n){
  elsa.cross <- data
  or.elsa.cross <- extract_RD(data = elsa.cross,
                              formula = f1,
                              start_var = disease_n) |>
    # filter(P_Value <0.05) |> 
    mutate(group = "Total",
           rank = rank(-RD, ties.method = "first")) |> 
    rbind(extract_RD(data = elsa.cross[elsa.cross$age.c == "<60",],
                     formula = f2,
                     start_var = disease_n) |> 
            #        filter(P_Value <0.05) |> 
            mutate(group = "<60",
                   rank = rank(-RD, ties.method = "first"))) |> 
    rbind(extract_RD(data = elsa.cross[elsa.cross$age.c == "60-74",],
                     formula = f2,
                     start_var = disease_n) |> 
            #       filter(P_Value <0.05) |> 
            mutate(group = "60-74",
                   rank = rank(-RD, ties.method = "first"))) |> 
    rbind(extract_RD(data = elsa.cross[elsa.cross$age.c == "75+",],
                     formula = f2,
                     start_var = disease_n) |> 
            #       filter(P_Value <0.05) |> 
            mutate(group = "75+",
                   rank = rank(-RD, ties.method = "first"))) |> 
    rbind(extract_RD(data = elsa.cross[elsa.cross$sex == "Male",],
                     formula = f3,
                     start_var = disease_n) |> 
            #       filter(P_Value <0.05) |> 
            mutate(group = "Male",
                   rank = rank(-RD, ties.method = "first"))) |> 
    rbind(extract_RD(data = elsa.cross[elsa.cross$sex == "Female",],
                     formula = f3,
                     start_var = disease_n) |> 
            #       filter(P_Value <0.05) |> 
            mutate(group = "Female",
                   rank = rank(-RD, ties.method = "first"))) 
    
}
