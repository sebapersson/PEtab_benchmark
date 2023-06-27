library(tidyverse)
library(ggthemes)
library(gt)
library(ggalluvial)


# General plotting parameters (plot using theme-tufte)
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
col_highlight = c("#D0C0B0", "#B6C2CC", "#BEAAB4", "#ECE9CD", "#0A3D6B", "#0D5C3D", "#812F02")
my_theme <- theme_hc(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                             plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
my_minimal <- theme_minimal(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                             plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0


read_data <- function(path_file)
{

  
  data_raw <- read_csv(path_file, col_types = cols()) 
  model_short = rep("", nrow(data_raw))
  for(i in 1:length(model_short)){
    iUse = str_locate(data_raw$model[i], "_[:alpha:]+_")
    model_short[i] = str_c(str_sub(data_raw$model[i], start=iUse[1]+1, end=iUse[2]-1), str_sub(data_raw$model[i], start=-4))
  }
  data_raw$model_short = model_short

  return(data_raw)
}


calc_winner_category <- function(data, tol=1e-6, sq_res=F)
{
  
  data_winner <- tibble()
  model_list <- unique(data$model)
  solver_type_list <- unique(data$solver_type)
  data <- data %>% filter(abstol == tol)
  
  for(i in 1:length(model_list)){
    for(j in 1:length(solver_type_list)){
      
      data_mod_i <- data %>%
        filter(model == model_list[i]) %>%
        filter(solver_type == solver_type_list[j]) %>%
        mutate(solver = as.factor(solver)) %>%
        group_by(solver) %>%
        summarise(median_time = median(runTime), 
                  median_sq = median(sqDiff), 
                  solver_type = median(solver_type), 
                  solver_lib = median(solver_lib), 
                  n_param = median(n_param), 
                  n_states = median(n_states), 
                  model_short = median(model_short), 
                  solverCategory = median(solverCategory))
      if(sq_res == FALSE){
        data_mod_i <- data_mod_i[which.min(data_mod_i$median_time)[1], ]
      }else{
        data_mod_i <- data_mod_i[which.min(data_mod_i$median_sq)[1], ]
      }
      
      data_winner <- data_winner %>%
        bind_rows(data_mod_i)
    }
  }
  
  return(data_winner)
}


plot_rank <- function(data, tol=1e-6, solver_type="all", sq_err=F)
{

  if(solver_type == "stiff"){
    data <- data %>%
      filter(solverType == "stiff")
    name_save = "Position_stiff"
  }else if(solver_type == "nonstiff"){
    data <- data %>%
      filter(solver_type != "stiff")
    name_save = "Position_not_stiff"
  }else{
    data <- data
    name_save = "Position"
  }
  
  data_winner = calc_winner(data, tol)
  
  if(sq_err == T){
    data_winner <- data_winner |> 
      rename("rank_plot" = "rank_sq") 
    ylab = "Ranking squared error"
    name_save1 <- str_c(name_save, "_sq")
    name_save2 <- str_c(name_save, "category_sq")
    dir_save <- str_c("Results/ODE_solvers/", as.character(tol), "/Sq_res/")
  }else if(sq_err == "Comb"){
    data_winner <- data_winner |> 
      mutate(rank_plot = (rank_sq + rank_time)*0.5)
    ylab = "Ranking averaged time and error"
    name_save1 <- str_c(name_save, "_comb")
    name_save2 <- str_c(name_save, "category_sq")
    dir_save <- str_c("Results/ODE_solvers/", as.character(tol), "/Res_comb/")
  }else{
    data_winner <- data_winner |> 
      mutate(rank_plot = rank_time)
    ylab = "Ranking squared error"
    name_save1 <- str_c(name_save, "_time")
    name_save2 <- str_c(name_save, "category_sq")
    dir_save <- str_c("Results/ODE_solvers/", as.character(tol), "/Time/")
  }
  
  posx = unique(data_winner$model_short[order(data_winner$n_states)])
  posy = data_winner |> group_by(solver) |> summarise(rank = median(rank_time))
  posy = posy$solver[order(posy$rank)]
  p1 = ggplot(data_winner, aes(x=model_short, y=solver)) + 
    geom_tile(aes(fill = rank_plot), color="white") + 
    scale_fill_viridis_c(direction = -1) +
    scale_x_discrete(limits = posx) + 
    scale_y_discrete(limits = rev(posy)) + 
    labs(x = "Model (increasing size ->)", y = ylab, title = "Rank time") + 
    my_minimal +
    theme(plot.background = element_rect(fill = "white"), 
          axis.text.x = element_text(angle=90))
  
  p2 <- ggplot(data_winner, aes(model_short)) + 
    scale_x_discrete(limits = posx) + 
    facet_grid(rows=vars(solverCategory), scales = "free_y", space="free_y") +
    geom_tile(aes(y=solver, fill = rank_plot), color="white") + 
    scale_fill_viridis_c(direction = -1) +
    labs(x = "Model (increasing size ->)", y = ylab, title = "Rank time") + 
    my_minimal +
    theme(plot.background = element_rect(fill = "white"), 
          axis.text.x = element_text(angle=90), 
          strip.placement = "outside")
  
  print(dir_save)
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  ggsave(str_c(dir_save, name_save1, ".svg"), p1, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*3.0)
  ggsave(str_c(dir_save, name_save2, ".svg"), p2, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*3.0)
}


calc_fail <- function(data, tol=1e-8)
{

  model_list <- unique(data$model)
  solver_list <- unique(data$solver)
  data_use <- data |> filter(abstol == tol)
  
  data_fail = data_use |> 
    group_by(solver, model_short, solverType, solverLib, solverCategory) |> 
    summarise(nFail = sum(is.nan(runTime) / 3)) |> 
    filter(nFail != 0)
    
  return(data_fail)
}


plot_fail <- function(data, tol)
{
  
  data_fail <- calc_fail(data, tol) 
  data1 = data_fail |> 
    group_by(model_short, solverType, solverCategory) |> 
    summarise(nFail = as.integer(sum(nFail)))
  data2 = data_fail |> 
    group_by(solver, solverType, solverCategory) |> 
    summarise(nFail = as.integer(sum(nFail)))
  data3 <- data_fail |> 
    mutate(category_ = case_when(solverLib == "Sundials" & solverType == "stiff" ~ "Sundials", 
                                 solverLib == "OrdinaryDiffEq" & solverType == "stiff" ~ "OrdinaryDiffEq stiff", 
                                 solverLib == "OrdinaryDiffEq" & solverType == "composite" ~ "OrdinaryDiffEq composite", 
                                 solverLib == "OrdinaryDiffEq" & solverType == "nonstiff" ~ "OrdinaryDiffEq nonstiff")) |> 
    group_by(solver, solverType, solverCategory, category_) |> 
    summarise(nFail = as.integer(sum(nFail)))
  data4 = data_fail |> 
      filter(solverType != "hint") |> 
    mutate(category_ = case_when(solverLib == "Sundials" & solverType == "stiff" ~ "Sundials", 
                                 solverLib == "OrdinaryDiffEq" & solverType == "stiff" ~ "OrdinaryDiffEq stiff", 
                                 solverLib == "OrdinaryDiffEq" & solverType == "composite" ~ "OrdinaryDiffEq composite", 
                                 solverLib == "OrdinaryDiffEq" & solverType == "nonstiff" ~ "OrdinaryDiffEq nonstiff")) |> 
    group_by(model_short, category_) |> 
    summarise(nFail = as.integer(sum(nFail))) 
  
  
  data1$solverType = factor(data1$solverType, levels = rev(c("nonstiff", "stiff", "hint", "composite")))
  pos1 = data1 |> group_by(model_short) |> summarise(nFail=sum(nFail)) 
  pos1 = pos1$model_short[order(pos1$nFail)]
  p1 = ggplot(data1, aes(model_short, nFail, fill = solverType)) + 
    geom_col(position=position_dodge(preserve = 'single')) + 
    scale_fill_manual(values = cbPalette[-1], name = "Solver type") + 
    scale_x_discrete(limits=pos1) + 
    scale_y_continuous(breaks = seq(from=1, by=2, to=max(data1$nFail))) +
    coord_flip() + 
    labs(y = "Number of integration failures", x = "", title = "Number of integration failures per model") + 
    my_minimal +
    theme(legend.position = "bottom", 
          plot.background = element_rect(fill = "white"))
  p2 = ggplot(data1, aes(model_short, nFail, fill = solverCategory)) + 
    geom_col(position=position_dodge(preserve = 'single')) + 
    scale_fill_manual(values = cbPalette, name = "Solver type") + 
    scale_x_discrete(limits=pos1) + 
    scale_y_continuous(breaks = seq(from=1, by=2, to=max(data1$nFail))) +
    coord_flip() + 
    labs(y = "Number of integration failures", x = "", title = "Number of integration failures per model") + 
    my_minimal +
    theme(legend.position = "bottom", 
          plot.background = element_rect(fill = "white"))
  
  pos2 = unique(data2$solver[order(data2$nFail)])
  p3 = ggplot(data2, aes(solver, nFail, fill = solverType)) + 
    geom_bar(stat="identity", position="dodge") + 
    scale_fill_manual(values = cbPalette[-1], name = "Solver type") + 
    scale_x_discrete(limits = pos2) +
    scale_y_continuous(breaks = seq(from=1, by=2, to=max(data2$nFail))) +
    coord_flip() + 
    labs(y = "Number of integration failures", x = "", title = "Number of integration failures (total of 27 models)") + 
    my_minimal + 
    theme(legend.position = "bottom", 
          plot.background = element_rect(fill = "white"))
  p4 = ggplot(data2, aes(solver, nFail, fill = solverCategory)) + 
    geom_bar(stat="identity", position="dodge") + 
    scale_fill_manual(values = cbPalette, name = "Solver type") + 
    scale_x_discrete(limits = pos2) +
    scale_y_continuous(breaks = seq(from=1, by=2, to=max(data2$nFail))) +
    coord_flip() + 
    labs(y = "Number of integration failures", x = "", title = "Number of integration failures (total of 27 models)") + 
    my_minimal + 
    theme(legend.position = "bottom", 
          plot.background = element_rect(fill = "white"))
  
  data3 <- data3 |> 
    drop_na() |> 
    filter(solverType != "hint")
  pos3 = unique(data3$solver[order(data3$nFail)])
  p5 = ggplot(data3, aes(solver, nFail, fill = category_)) + 
    geom_bar(stat="identity", position="dodge") + 
    scale_fill_manual(values = cbPalette[-1], name = "Solver type") + 
    scale_x_discrete(limits = pos3) +
    scale_y_continuous(breaks = seq(from=1, by=2, to=max(data3$nFail))) +
    coord_flip() + 
    labs(y = "Number of integration failures", x = "", title = "Number of integration failures (total of 27 models)") + 
    my_minimal + 
    theme(legend.position = "bottom", 
          plot.background = element_rect(fill = "white"))
  
  p6 = ggplot(data4, aes(model_short, nFail, fill = category_)) + 
    geom_bar(stat="identity", position="dodge") + 
    scale_fill_manual(values = cbPalette[-1], name = "Solver type") + 
    scale_x_discrete(limits = pos1) +
    scale_y_continuous(breaks = seq(from=1, by=2, to=max(data3$nFail))) +
    coord_flip() + 
    labs(y = "Number of integration failures", x = "", title = "Number of integration failures (total of 27 models)") + 
    my_minimal + 
    theme(legend.position = "bottom", 
          plot.background = element_rect(fill = "white"))
  
  dir_save <- str_c("Results/ODE_solvers/", as.character(tol), "/")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  ggsave(str_c(dir_save, "Solvers_fail_solver.svg"), p1, width = BASE_WIDTH*2, height = BASE_HEIGHT*2.8)
  ggsave(str_c(dir_save, "Solvers_fail_solver_category.svg"), p2, width = BASE_WIDTH*2, height = BASE_HEIGHT*2.8)
  ggsave(str_c(dir_save, "Solvers_fail_model.svg"), p3, width = BASE_WIDTH*2, height = BASE_HEIGHT*2.8)
  ggsave(str_c(dir_save, "Solvers_fail_model_category.svg"), p4, width = BASE_WIDTH*2, height = BASE_HEIGHT*2.8)
  ggsave(str_c(dir_save, "Solvers_fail_model_lib.svg"), p5, width = BASE_WIDTH*2, height = BASE_HEIGHT*2.8)
  ggsave(str_c(dir_save, "Solvers_fail_model_x_lib.svg"), p6, width = BASE_WIDTH*2, height = BASE_HEIGHT*2.8)
}


calc_winner <- function(data, tol=1e-6, na_rm=T)
{
  
  return_first = function(value) return(value[1])
  
  data_winner <- tibble()
  data_use <- data |>  
    filter(abstol == tol) 
  model_list <- unique(data$model_short)
  n_solvers = length(unique(data$solver))
  
  for(i in 1:length(model_list)){
    
    data_model_i <- data_use |> 
      filter(model_short == model_list[i]) |> 
      group_by(solver, model_short) |> 
      summarise(median_time = median(runTime, na.rm = na_rm), 
                median_sq = median(sqDiff, na.rm = na_rm), 
                solver_type = return_first(solverType), 
                solver_lib = return_first(solverLib), 
                n_param = median(nParam), 
                n_states = median(nStates), 
                solverCategory = return_first(solverCategory))
    # Add ranked index on time or sq-error
    data_model_i <- data_model_i[order(data_model_i$median_time), ]
    data_model_i$rank_time <- 1:nrow(data_model_i)
    data_model_i <- data_model_i[order(data_model_i$median_sq), ]
    data_model_i$rank_sq <- 1:nrow(data_model_i)
    
    # Account for integration failures     
    data_model_i$rank_time[is.na(data_model_i$median_time)] = NaN
    data_model_i$rank_sq[is.na(data_model_i$median_sq)] = NaN
    
    data_winner <- bind_rows(data_winner, data_model_i)
  }
  
  return(data_winner)
}


plot_winner <- function(data, tol, sq_res=F){
  
  data_winner <- calc_winner(data, tol=tol, na_rm=F)
  
  if(sq_res == T){
    data_winner <- data_winner %>%
      rename("crit_plot" = "median_sq")
    ylab = "Squared error"
  }else{
    data_winner <- data_winner %>%
      rename("crit_plot" = "median_time")
    ylab = "Time [s]"
  }
  
  # Line plot for stiff solvers 
  data_stiff <- data_winner |> 
    filter(solver_type == "stiff")
  data_stiff <- data_stiff |> mutate(crit_plot = case_when(is.infinite(crit_plot) ~ max(crit_plot, na.rm = T)+0.5, 
                                                           is.nan(crit_plot) ~ max(crit_plot, na.rm = T)+0.5, 
                                                           is.na(crit_plot) ~ max(crit_plot, na.rm = T)+0.5, 
                                                           T ~ crit_plot))
  
  pos <- unique(data_stiff$model_short[order(data_stiff$n_states)])
  stiff_solvers_plot <- c("Rodas5P", "CVODE_BDF_Dense", "RadauIIA3", "QNDF")
  solvers_color <- data_stiff |> filter(solver %in% stiff_solvers_plot)
  p1 <- ggplot(data_stiff, aes(x = model_short, crit_plot)) + 
    geom_point(aes(group=solver), stat="summary", fun=median, color = cbPalette[1]) + 
    stat_summary(aes(group=solver), fun=median, geom="line", color = cbPalette[1]) +
    geom_point(data=solvers_color, mapping=aes(color=solver, group=solver), stat="summary", fun=median, size=3.0) + 
    geom_line(data=solvers_color, mapping=aes(color=solver, group=solver), stat="summary", fun=median, linewidth=2.0) + 
    scale_y_log10() + 
    labs(x = "", y = ylab, title = "Stiff solvers") +
    scale_color_manual(values = cbPalette[c(5, 6, 7, 8)], name = "Selected solvers") + 
    scale_x_discrete(limits=pos) + 
    my_minimal +
    theme(axis.text.x = element_text(angle = 90), 
          legend.position = "bottom", 
          plot.background = element_rect(fill = "white"))
  
  data_non_stiff <-  data_winner %>%
    filter(solver_type != "stiff") |> 
    filter(solver_type != "hint") |> 
    mutate(crit_plot = case_when(is.nan(crit_plot) ~ max(crit_plot, na.rm = T)+0.5, 
                                                           is.nan(crit_plot) ~ max(crit_plot, na.rm = T)+0.5, 
                                                           is.na(crit_plot) ~ max(crit_plot, na.rm = T)+0.5, 
                                                           T ~ crit_plot))
  nonstiff_solvers_plot <- c("Vern7", "BS3", "Vern7Rodas5", "Tsit5", "QNDF")
  solvers_color <- data_non_stiff |> filter(solver %in% nonstiff_solvers_plot)
  p2 <- ggplot(data_non_stiff, aes(x = model_short, crit_plot)) + 
    geom_point(aes(group=solver), stat="summary", fun=median, color = cbPalette[1]) + 
    stat_summary(aes(group=solver), fun=median, geom="line", color = cbPalette[1]) +
    geom_point(data=solvers_color, mapping=aes(color=solver, group=solver), stat="summary", fun=median, size=3.0) + 
    geom_line(data=solvers_color, mapping=aes(color=solver, group=solver), stat="summary", fun=median, linewidth=2.0) + 
    scale_y_log10() + 
    labs(x = "", y = ylab, title = "Stiff solvers") +
    scale_color_manual(values = cbPalette[-1], name = "Selected solvers") + 
    scale_x_discrete(limits=pos) + 
    my_minimal +
    theme(axis.text.x = element_text(angle = 90), 
          legend.position = "bottom", 
          plot.background = element_rect(fill = "white"))
  
  # Library and solver type. Extract for each model the best performer 
  data_lib_ <- data_winner |> 
    filter(solver_type != "hint") |> 
    mutate(category_ = case_when(solver_lib == "Sundials" & solver_type == "stiff" ~ "Sundials", 
                                 solver_lib == "OrdinaryDiffEq" & solver_type == "stiff" ~ "OrdinaryDiffEq stiff", 
                                 solver_lib == "OrdinaryDiffEq" & solver_type == "composite" ~ "OrdinaryDiffEq composite", 
                                 solver_lib == "OrdinaryDiffEq" & solver_type == "nonstiff" ~ "OrdinaryDiffEq nonstiff")) 
  data_lib <- tibble()
  for(i in 1:length(unique(data_lib_$category_))){
    for(j in 1:length(unique(data_winner$model_short))){
      data_tmp <- data_lib_ |> 
        filter(category_ == unique(data_lib_$category_)[i]) |> 
        filter(model_short == unique(data_winner$model_short)[j])
      if(sq_res == TRUE){
        if(is.infinite(min(data_tmp$rank_sq, na.rm = T))){
          data_tmp = data_tmp[1, ]
        }else{
          data_tmp <- data_tmp |> filter(rank_sq == min(rank_sq, na.rm = T))
        }
      }else{
        if(is.infinite(min(data_tmp$rank_time, na.rm = T))){
          data_tmp = data_tmp[1, ]
        }else{
        data_tmp <- data_tmp |> filter(rank_time == min(rank_time, na.rm = T))
        }
      }
      data_lib <- bind_rows(data_lib, data_tmp)
    }
  }
  data_lib <- data_lib |> 
    mutate(crit_plot = case_when(is.infinite(crit_plot) ~ max(crit_plot, na.rm = T)+0.9, 
                                 is.nan(crit_plot) ~ max(crit_plot, na.rm = T)+0.9, 
                                 is.na(crit_plot) ~ max(crit_plot, na.rm = T)+0.9, 
                                 T ~ crit_plot))
  
  p3 <- ggplot(data_lib, aes(x = model_short, crit_plot, group=category_, color = category_)) + 
    geom_point(stat="summary", fun=median, size = 3.0) +
    stat_summary(fun=median, geom="line", size = 1.5) +
    scale_y_log10() + 
    labs(x = "", y = ylab, title = "Performance different solver categories") +
    scale_color_manual(values = cbPalette[-1]) + 
    scale_x_discrete(limits=pos) +  
    my_minimal +    
    theme(axis.text.x = element_text(angle = 90), 
          legend.position = "bottom", 
          plot.background = element_rect(fill = "white"))
  
  # For avg distance 
  if(sq_res == F){
    model_small1 <- data_lib |> filter(n_states <= 16 & category_ == "Sundials")
    model_small2 <- data_lib |> filter(n_states <= 16 & category_ == "OrdinaryDiffEq stiff")
    avg_time1 <- mean(model_small1$crit_plot / model_small2$crit_plot) 
    model_large1 <- data_lib |> filter(n_states > 16 & category_ == "Sundials")
    model_large2 <- data_lib |> filter(n_states > 16 & category_ == "OrdinaryDiffEq stiff")
    avg_time2 <- mean(model_large1$crit_plot / model_large2$crit_plot) 
    print(sprintf("Avg. speedup small models = %.3f", avg_time1))
  }
  
  # Build a table with results per model 
  my_tab <- data_lib |> 
    select(model_short, crit_plot, solver, category_) |> 
    filter(crit_plot != max(crit_plot)) |> 
    gt(rowname_col = "model_short") 
  for(model in unique(data_lib$model_short)){
    my_tab <- my_tab |> 
      tab_row_group(label = model, rows = starts_with(model))
  }
  my_tab <- my_tab |> 
    cols_label(solver= "ODE solver", crit_plot=ylab, category_ = "Solver cateogry") |> 
    fmt(columns = starts_with("crit_plot"), fns = function(x) sprintf("%.2e", x)) |> 
    fmt(columns = starts_with("model_short"), fns = function(x) sprintf("")) 
  
  if(sq_res == T){
    dir_save <- str_c("Results/ODE_solvers/", as.character(tol), "/Sq_res/")
  }else{
    dir_save <- str_c("Results/ODE_solvers/", as.character(tol), "/Time/")
  }
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  
  ggsave(str_c(dir_save, "Stiff_solvers.svg"), p1, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Non_stiff_solvers.svg"), p2, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Solver_category.svg"), p3, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  gtsave(my_tab, str_c(dir_save, "Table_result.html"))
}


plot_linear_solvers <- function(data, tol)
{

  solvers_test <- c("QNDF", "RadauIIA5", "Rodas5", "Rodas4P", "CVODE_BDF")
  
  for(i in 1:length(solvers_test)){
    data_solver <- data |> 
      filter(abstol == tol) |> 
      filter(grepl(str_c("^",solvers_test[i]), solver)) |> 
      group_by(solver, model_short) |> 
      summarise(run_time = median(runTime), 
                sd_time = sd(runTime)) |> 
      mutate(linear_solver = case_when(str_sub(solver, -5, -1) == "KLU_S" ~ "KLU Sparse", 
                                       str_sub(solver, -5, -1) != "KLU_S" & str_sub(solver, -7, -1) != "GMRES_S" &  str_sub(solver, -2, -1) == "_S" ~ "Sparse", 
                                       str_sub(solver, -7, -1) == "GMRES_S" ~ "GMRES Sparse", 
                                       str_sub(solver, -5, -1) == "GMRES" ~ "GMRES", 
                                       str_sub(solver, -6, -1) == "FastLU" ~ "FastLU", 
                                       str_sub(solver, -5, -1) == "RFLUF" ~ "RFLUF", 
                                       solver == "CVODE_BDF_default" ~ "Standard", 
                                       str_sub(solver, -11, -1) == "LapackDense" ~ "Lapack", 
                                       str_sub(solver, -5, -1) == "Dense" ~ "Dense", 
                                       T ~ "Standrad")) |> 
      filter(linear_solver != "GMRES") |> 
      filter(solver != "Rodas5P")
    
    pos <- unique(data$model_short[order(data$nStates)])
    p <- ggplot(data_solver, aes(x = model_short, run_time, color=linear_solver)) + 
      geom_line(aes(group=solver), stat="summary", fun=median, linewidth=2.0) + 
      scale_y_log10() + 
      labs(x = "", y = "Run time [s]", title = str_c(solvers_test[i], " : Different linear-solvers")) +
      scale_color_manual(values = cbPalette[-1], name = "Selected linear solvers") + 
      scale_x_discrete(limits=pos) + 
      my_minimal +
      theme(axis.text.x = element_text(angle = 90), 
            legend.position = "bottom", 
            plot.background = element_rect(fill = "white"))
    
    dir_save <- str_c("Results/ODE_solvers/", as.character(tol), "/Time/Linear_solvers/")
    if(!dir.exists(dir_save)) dir.create(dir_save)
    ggsave(str_c(dir_save, solvers_test[i], ".svg"), p, width = BASE_WIDTH*2.0, height = BASE_HEIGHT*2.0)
  }
}


data1 <- read_data("Master-Thesis/Intermediate/Benchmarks/ODE_solvers/All_models.csv") |> 
  mutate(solverType =  case_when(solverType == "nonStiff" ~ "nonstiff", 
                                 T ~ solverType)) 
data2 <- read_data("Master-Thesis/Intermediate/Benchmarks/ODE_solvers/All_models_sparse_jacobian.csv") |> 
  mutate(solverType =  case_when(solverType == "nonStiff" ~ "nonstiff", 
                                 T ~ solverType)) 
data <- bind_rows(data1, data2) |> 
  filter(solver != "RadauIIA5_GMRES_S") |> 
  mutate(linear_solver = case_when(str_sub(solver, -5, -1) == "KLU_S" ~ "KLU_Sparse", 
                                   str_sub(solver, -5, -1) != "KLU_S" & str_sub(solver, -2, -1) == "_S" ~ "Sparse", 
                                   str_sub(solver, -7, -1) == "GMRES_S" ~ "Sparse GMRES", 
                                   str_sub(solver, -6, -1) == "FastLU" ~ "FastLU", 
                                   str_sub(solver, -5, -1) == "RFLUF" ~ "RFLUF", 
                                   str_sub(solver, -6, -1) == "_GMRES" ~ "Sparse GMRES",
                                   str_sub(solver, -6, -1) == "_Dense" ~ "Dense",
                                   str_sub(solver, -11, -1) == "LapackDense" ~ "Lapack Dense",
                                   T ~ "Standard")) |> 
  filter(linear_solver == "Standard")

tol_list <- c(1e-6, 1e-8, 1e-16)
for(tol in tol_list){
  plot_fail(data, tol)
  plot_rank(data, tol, solver_type = "all", sq_err = F)
  plot_rank(data, tol, solver_type = "stiff", sq_err = F)
  plot_rank(data, tol, solver_type = "nonstiff", sq_err = F)
  plot_rank(data, tol, solver_type = "all", sq_err = T)
  plot_rank(data, tol, solver_type = "stiff", sq_err = T)
  plot_rank(data, tol, solver_type = "nonstiff", sq_err = T)
  plot_winner(data, tol, sq_res = F)
  plot_winner(data, tol, sq_res = T)
  plot_linear_solvers(data, tol)
}


# -----------------------------------------------------------------------------------------------------------------
# Large models 
# -----------------------------------------------------------------------------------------------------------------
data_large <- read_data("Master-Thesis/Intermediate/Benchmarks/ODE_solvers/Large_models_random_p.csv") |> 
  mutate(solverType =  case_when(solverType == "nonStiff" ~ "nonstiff", 
                                 T ~ solverType)) |> 
  group_by(model, solver, solverLib, model_short, index_param) |> 
  summarise(runTime = median(runTime)) |> 
  mutate(is_stiff = str_sub(solver, -2, -1) == "_S") |> 
  mutate(model_short = case_when(is.na(model_short) ~ "Smith2013", 
                                 T ~ model_short)) |> 
  filter(str_sub(solver, 1, 5) != "Rodas") |> 
  filter(str_sub(solver, 1, 5) != "TRBDF") |> 
  filter(str_sub(solver, 1, 9) != "KenCarp47") |> 
  mutate(linear_solver = case_when(str_sub(solver, -5, -1) == "KLU_S" ~ "KLU_Sparse", 
                                   str_sub(solver, -5, -1) != "KLU_S" & str_sub(solver, -2, -1) == "_S" ~ "Sparse", 
                                   str_sub(solver, -6, -1) == "FastLU" ~ "FastLU", 
                                   str_sub(solver, -5, -1) == "RFLUF" ~ "RFLUF", 
                                   T ~ "Standrad")) |> 
  mutate(solver_short = case_when(str_sub(solver, -5, -1) == "KLU_S" ~ str_sub(solver, 1, -7), 
                                  str_sub(solver, -5, -1) != "KLU_S" & str_sub(solver, -2, -1) == "_S" ~ str_sub(solver, 1, -3), 
                                  str_sub(solver, -6, -1) == "FastLU" ~ str_sub(solver, 1, -8), 
                                  str_sub(solver, -5, -1) == "RFLUF" ~ str_sub(solver, 1, -7), 
                                  solver == "CVODE_BDF_default" ~ "CVODE_BDF",
                                  T ~ solver))

data_large$model_short <- factor(data_large$model_short, levels = c("Bachmann2011", "Isensee2018", "Lucarelli2018", "Smith2013", "Chen2009"))
data_large$linear_solver <- factor(data_large$linear_solver, levels = rev(c("Standrad", "RFLUF", "FastLU", "Sparse", "KLU_Sparse")))

## First we look at the nominal (reported) values
data_nominal <- data_large |> 
  filter(index_param == 1) 
data_best <- data_nominal |> 
  group_by(model) |> 
  summarise(min_time = min(runTime, na.rm = T))
data_nominal <- inner_join(data_nominal, data_best, by = c("model"))

ggplot(data_nominal, aes(solver_short, runTime / min_time, fill = linear_solver)) + 
  geom_bar(stat="identity", position = position_dodge2(preserve = "single", padding = 0)) + 
  scale_y_continuous(expand=c(0, 0, 0.1, 0.1), breaks = seq(from=1, by = 2, to = 20)) +
  facet_wrap(~model_short, scale="free_x") + 
  coord_flip() +
  geom_hline(yintercept = 1.0) +
  scale_fill_manual(values = cbPalette[-1], name = "Linear solver") + 
  labs(y = "Run time normalised by the fastest solver for each model", x = "") + 
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")

# Look wider 
data_large_plot <- data_large |> 
  group_by(solver_short, linear_solver, model_short) |> 
  summarise(median_time = median(runTime, na.rm = T),
            min_time = min(runTime, na.rm = T), 
            max_time = max(runTime, na.rm = T),
            n_fail = sum(is.na(runTime), na.rm = T))

data_best <- data_large_plot |> 
  group_by(model_short) |> 
  summarise(min_median_time = min(median_time, na.rm = T))

data_large_plot <- inner_join(data_large_plot, data_best, by = c("model_short")) |> 
  mutate(time_norm = median_time / min_median_time )

p_large_norm = ggplot(data_large_plot, aes(solver_short, time_norm, fill = linear_solver)) + 
  geom_bar(stat="identity", position = position_dodge2(preserve = "single", padding = 0)) + 
  geom_text(aes(x=solver_short, y = time_norm * 1.06, label = sprintf("%d", n_fail)), 
            position = position_dodge2(width=1.0)) +
  facet_wrap(~model_short, scale="free_x") + 
  coord_flip() +
  scale_y_continuous(expand=c(0, 0, 0.1, 0.1), breaks = seq(from=1, by = 2, to = 30)) +
  scale_fill_manual(values = cbPalette[-1], name = "Linear solver") + 
  labs(y = "Run time normalised by the fastest solver for each model", x = "") + 
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")

data_fail <- data_large |> 
  group_by(solver_short, linear_solver, model_short) |> 
  summarise(max_time = max(runTime, na.rm = T),
            n_fail = sum(is.na(runTime), na.rm = T))

p_large = ggplot(data_large, aes(solver_short, runTime, fill = linear_solver)) + 
  geom_boxplot(position = position_dodge2(preserve = "single")) + 
  geom_text(data=data_fail, mapping=aes(x=solver_short, y = max_time * 1.02, label = sprintf("%d", n_fail)), 
            position = position_dodge2(width=1.0)) +
  facet_wrap(~model_short, scale="free_x") + 
  coord_flip() +
  scale_fill_manual(values = cbPalette[-1], name = "Linear solver") + 
  scale_y_log10() + 
  labs(y = "Run time across 100 random parameter vectors [s]", x = "") + 
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")

ggplot(data_fail, aes(solver_short, n_fail, fill = linear_solver)) + 
  geom_bar(stat="identity")

data_tmp = data_large |> 
  mutate(fail_cate = case_when(is.na(runTime) ~ "Yes", 
                               T ~ "No")) |> 
  mutate(fail_cate = factor(fail_cate, levels = c("No", "Yes")))
data_tmp1 = data_tmp |>      
    filter(fail_cate == "Yes")
data_tmp2 = data_tmp |>      
    filter(fail_cate == "No")
data_tmp = bind_rows(data_tmp1, data_tmp2)
         
p = ggplot(data_tmp, aes(axis1=solver_short, 
                         axis2=model_short)) +
  geom_alluvium(aes(fill=fail_cate), segments=1000, width = 1/12, alpha=1.0) +
  geom_stratum(width = 1/12, fill = cbPalette[1], color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = col_highlight[c(1, 5)]) +
  theme_void() + 
  scale_x_discrete(expand=c(0, 0)) + 
  theme(legend.position = "none")


p <- ggplot(data_fail, aes(solver_short, n_fail, fill = linear_solver)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~model_short) +
  scale_fill_manual(values = cbPalette[-1], name = "Linear solver") + 
  scale_y_continuous(expand=c(0, 0)) + 
  theme_bw(base_size=22) + 
  labs(x = "", y = "Number of integration failures", 
       title = "Number of integration failures", 
       subtitle = "50 random parameter vectors") + 
  coord_flip() + 
  theme(legend.position = "bottom")

  
ggsave("Rplot.svg", p, dpi=900, width = BASE_WIDTH, height = BASE_HEIGHT)

ggsave("Results/ODE_solvers/Integration_failures.png", dpi=300, width = BASE_WIDTH*2.0, height = BASE_HEIGHT*2.0)
ggsave("Results/ODE_solvers/Large_norm.svg", p_large_norm, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*2.0)
ggsave("Results/ODE_solvers/Large.svg", p_large, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*2.0)

# -----------------------------------------------------------------------------------------------------------------
# Random parameters 
# -----------------------------------------------------------------------------------------------------------------
dir_data = "Master-Thesis/Intermediate/Benchmarks/ODE_solvers/"
data_random_parameters = read_data(str_c(dir_data, "Random_parameters.csv")) |> 
  mutate(solverType = case_when(solverType == "nonStiff" ~ "nonstiff", 
                                T ~ solverType))

data_plot = data_random_parameters |> 
  filter(abstol == 1e-8) |> 
  mutate(solver = case_when(solver == "CVODE_BDF_default" ~ "CVODE_BDF", 
                            T ~ solver))
data_plot = data_plot |> 
  mutate(label = case_when(model_short %in% c("Perelson1996", "Zhao2020", "Crauste2017", "Fiedler2016", "Bruno2016", "Okuonghae2020") ~ "Non-stiff solver fastest", 
                           model_short %in% c("Schwen2014", "Bachmann2011", "Lucarelli2018") ~ "Stiff solver fastest", 
                           T ~ "Has preequlibrium")) |> 
  filter(model != "model_Brannmark_JBC2010") |> 
  mutate(category_ = case_when(solverLib == "Sundials" & solverType == "stiff" ~ "Sundials", 
                               solverLib == "OrdinaryDiffEq" & solverType == "stiff" ~ "OrdinaryDiffEq stiff", 
                               solverLib == "OrdinaryDiffEq" & solverType == "composite" ~ "OrdinaryDiffEq composite", 
                               solverLib == "OrdinaryDiffEq" & solverType == "nonstiff" ~ "OrdinaryDiffEq nonstiff")) |> 
  filter(solver != "Rodas4P") |> 
  filter(solver != "Vern6") |> 
  filter(model_short %in% c("Perelson1996", "Zhao2020", "Crauste2017", "Fiedler2016", "Bruno2016", 
                            "Okuonghae2020", "Schwen2014", "Bachmann2011", "Bertozzi2020", 
                            "Giordano2020", "Borghans1997", "Weber2015"))


data_main <- data_plot |> 
  filter(model_short %in% c("Perelson1996", "Borghans1997", "Crauste2017", "Bruno2016", 
                            "Bachmann2011", "Weber2015", "Schwen2014"))
levels_use <- data_main |> 
  group_by(model_short) |> 
  summarise(n_states = mean(nStates))
levels_use_ <- levels_use$model_short[order(levels_use$n_states)]
data_main$model_short = factor(data_main$model_short, levels = levels_use_)

data_not_main <- data_plot |> 
  filter(!model_short %in% c("Perelson1996", "Borghans1997", "Crauste2017", "Bruno2016", 
                             "Bachmann2011", "Weber2015", "Schwen2014"))
levels_use <- data_not_main |> 
  group_by(model_short) |> 
  summarise(n_states = mean(nStates))
levels_use_ <- levels_use$model_short[order(levels_use$n_states)]
data_not_main$model_short = factor(data_not_main$model_short, levels = levels_use_)
data_fail_not = data_not_main |> 
  group_by(solver, model_short, category_, label) |> 
  summarise(n_fail = sum(is.nan(runTime)), 
            runTime = max(runTime, na.rm = T) + max(runTime, na.rm = T)*0.2) |> 
  mutate(runTime = case_when(is.nan(runTime) ~ 0.01, 
                              T ~ runTime))


data_fail = data_main |> 
  group_by(solver, model_short, category_, label) |> 
  summarise(n_fail = sum(is.nan(runTime)), 
            runTime = max(runTime, na.rm = T) + max(runTime, na.rm = T)*0.2) |> 
  mutate(runTime = case_when(is.nan(runTime) ~ 0.01, 
                              T ~ runTime))

p1 = ggplot(data_main, aes(solver, runTime)) + 
  geom_violin(aes(fill = category_), draw_quantiles = 0.5, linewidth=1.0) +
  geom_jitter(width = 0.1, size=0.2) +
  geom_text(data=data_fail, mapping=aes(x=solver, y=runTime, label=sprintf("%d", n_fail)), size=8.0) + 
  facet_wrap(~model_short, scales="free_x", nrow=4, ncol=2) + 
  scale_y_log10() +
  labs(x = "", y = "Run time [s]", title = "Run time for 100 random paramter vectors", 
       subtitle = "On average stiff solvers have less variabillity") + 
  scale_fill_manual(values = cbPalette[-1], name = "Solver type") +
  theme_bw(base_size = 20) +
  coord_flip() + 
  theme(legend.position = "bottom")


p_sir <- ggplot(data_not_main, aes(solver, runTime)) + 
  geom_violin(aes(fill = category_), draw_quantiles = 0.5, linewidth=1.0) +
  geom_jitter(width = 0.1, size=0.2) +
  geom_text(data=data_fail_not, mapping=aes(x=solver, y=runTime, label=sprintf("%d", n_fail)), size=8.0) + 
  facet_wrap(~model_short, scales="free_x", nrow=1, ncol=4) + 
  scale_y_log10() +
  labs(x = "", y = "Run time [s]", title = "Run time for 100 random paramter vectors", 
       subtitle = "On average stiff solvers have less variabillity") + 
  scale_fill_manual(values = cbPalette[-1], name = "Solver type") +
  theme_bw(base_size = 20) +
  coord_flip() + 
  theme(legend.position = "bottom")


p2 <- ggplot(data_fail, aes(solver, n_fail, fill = category_)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~model_short) + 
  labs(x = "", y = "Percentage integration failures", title = "Percentage integration failure for 100 random parameter vectors", 
       subtitle = "On average nonstiff solvers have more integration failures") + 
  scale_fill_manual(values = cbPalette[-1], name = "Solver type") +
  theme_bw(base_size = 16) + 
  scale_y_continuous(expand=expansion(add = c(0, 5)), limits = c(0, 100)) +
  geom_text(x = "Vern7Rodas4P", y = 70, aes(label=label)) + 
  coord_flip() + 
  theme(legend.position = "bottom")


ggsave("Results/ODE_solvers/Random_parameter_run_time.svg", p1, width = BASE_WIDTH*2.0, height = BASE_HEIGHT*4.0)
ggsave("Results/ODE_solvers/Random_parameter_run_time_sir.svg", p_sir, width = BASE_WIDTH*3.0, height = BASE_HEIGHT*1.0)
ggsave("Results/ODE_solvers/Random_parameter_fail.png", p2, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*3.5, dpi=300)
