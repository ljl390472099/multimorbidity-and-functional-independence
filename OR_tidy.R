OR_tidy <- function(data, f1, f2, disease_n, sig = F){
  # sig设置是否需要按显著性筛选，F表示不筛选
  elsa.cross <- data
  if (sig == F){
    or.elsa.cross <- extract_OR(data = elsa.cross,
                                formula = f1,
                                extract_n = disease_n) |>
      # filter(P_Value <0.05) |> 
      mutate(group = "Total",
             rank = rank(-OR, ties.method = "first")) |> 
      rbind(extract_OR(data = elsa.cross[elsa.cross$age.c == "<60",],
                       formula = f2,
                       extract_n = disease_n) |> 
              #        filter(P_Value <0.05) |> 
              mutate(group = "<60",
                     rank = rank(-OR, ties.method = "first"))) |> 
      rbind(extract_OR(data = elsa.cross[elsa.cross$age.c == "60-74",],
                       formula = f2,
                       extract_n = disease_n) |> 
              #       filter(P_Value <0.05) |> 
              mutate(group = "60-74",
                     rank = rank(-OR, ties.method = "first"))) |> 
      rbind(extract_OR(data = elsa.cross[elsa.cross$age.c == "75+",],
                       formula = f2,
                       extract_n = disease_n) |> 
              #       filter(P_Value <0.05) |> 
              mutate(group = "75+",
                     rank = rank(-OR, ties.method = "first")))
    
  } else{
    or.elsa.cross <- extract_OR(data = elsa.cross,
                                formula = f1,
                                extract_n = disease_n) |>
       filter(P_Value <0.05) |> 
      mutate(group = "Total",
             rank = rank(-OR, ties.method = "first")) |> 
      rbind(extract_OR(data = elsa.cross[elsa.cross$age.c == "<60",],
                       formula = f2,
                       extract_n = disease_n) |> 
                       filter(P_Value <0.05) |> 
              mutate(group = "<60",
                     rank = rank(-OR, ties.method = "first"))) |> 
      rbind(extract_OR(data = elsa.cross[elsa.cross$age.c == "60-74",],
                       formula = f2,
                       extract_n = disease_n) |> 
                       filter(P_Value <0.05) |> 
              mutate(group = "60-74",
                     rank = rank(-OR, ties.method = "first"))) |> 
      rbind(extract_OR(data = elsa.cross[elsa.cross$age.c == "75+",],
                       formula = f2,
                       extract_n = disease_n) |> 
                       filter(P_Value <0.05) |> 
              mutate(group = "75+",
                     rank = rank(-OR, ties.method = "first")))
    

  }
  return(or.elsa.cross)
}
