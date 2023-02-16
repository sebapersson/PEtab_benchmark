library(tidyverse)
library(ggthemes)
library(gt)
library(rhdf5)


# General plotting parameters (plot using theme-tufte)
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
col_highlight = c("#D0C0B0", "#B6C2CC", "#BEAAB4", "#ECE9CD", "#0A3D6B", "#0D5C3D", "#812F02")
my_theme <- theme_hc(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                             plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
my_minimal <- theme_minimal(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                                    plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0


get_niter <- function(stat_file, stat_names)
{
  n_iter <- rep(0, length(stat_names))
  for(j in 1:length(stat_names)){
    n_iter[j] <- length(stat_file[[stat_names[j]]]$fval)
  }
  return(n_iter)
}


get_name_order <- function(stat_file, stat_names)
{
  fval = rep(0, length(stat_names))
  for(i in 1:length(fval)){
    fval_tmp <- stat_file[[stat_names[i]]]$fval
    fval[i] <- fval_tmp[length(fval_tmp)]
  }
  return(order(fval))
}


file <- "/home/sebpe/Dropbox/PhD/Projects/PEtab_benchmark/pypesto_benchmark/results/Fujita_SciSignal2010__fides.hessian=FIM__1000.hdf5"
pyPestoResFile <- h5read(file, "/")
pyPestoResNames <- names(pyPestoResFile)


get_pyPesto_results <- function(model_name)
{
  
  pyPestoFiles <- list.files("pypesto_benchmark/results")[str_detect(list.files("pypesto_benchmark/results"), str_c(model_name, ".+", ".csv"))]
  pyPestoStatFiles <- str_c("pypesto_benchmark/stats/",  list.files("pypesto_benchmark/stats")[str_detect(list.files("pypesto_benchmark/stats"), str_c(model_name, ".+", ".hdf5"))])
  pyPestoData <- tibble()
  
  for(i in 1:length(pyPestoFiles)){
  
    pyPestoData_ <- read_csv(str_c("pypesto_benchmark/results/", pyPestoFiles[i]), col_names = F, col_types = cols()) |> 
      filter(!is.infinite(X1))
    names(pyPestoData_) <- c("Cost", "Run_time")
    pyPestoStatFile <- h5read(pyPestoStatFiles[i], "/")
    pyPestoStatsNames <- names(pyPestoStatFile)
    pyPestoStatsNames <- pyPestoStatsNames[get_name_order(pyPestoStatFile, pyPestoStatsNames)]
    
    if(str_detect(pyPestoFiles[i], "fides.hessian=FIM") == T){
      pyPestoData_ <- pyPestoData_ |> 
        mutate(Alg = "PyPesto Fides FIM", 
               Solver = "AMICI", 
               N_iter = get_niter(pyPestoStatFile, pyPestoStatsNames))
    }else if(str_detect(pyPestoFiles[i], "fides.hessian=BFGS") == T){
      pyPestoData_ <- pyPestoData_ |> 
        mutate(Alg = "PyPesto Fides BFGS", 
               Solver = "AMICI", 
               N_iter = get_niter(pyPestoStatFile, pyPestoStatsNames))
    }else if(str_detect(pyPestoFiles[i], "fides__") == T){
      pyPestoData_ <- pyPestoData_ |> 
        mutate(Alg = "PyPesto Fides", 
               Solver = "AMICI", 
               N_iter = get_niter(pyPestoStatFile, pyPestoStatsNames))
    }else{
      pyPestoData_ <- tibble()
    }
    h5closeAll()
    pyPestoData <- bind_rows(pyPestoData, pyPestoData_)
  }
  
  return(pyPestoData)
}


get_rank <- function(data)
{
  data_ret <- tibble()
  all_solvers <- unique(data$Alg)
  for(i in 1:length(all_solvers)){
    data_solver <- data |> 
      filter(Alg == all_solvers[i]) 
    data_solver <- data_solver[order(data_solver$Cost), ] |> 
      mutate(rank = 1:nrow(data_solver))
    data_ret <- bind_rows(data_ret, data_solver)
  }
  return(data_ret)
}


process_model <- function(model_name, ylim=NULL, jl_solver=NULL)
{
  
  dir_julia <- str_c("Master-Thesis/Intermediate/Benchmarks/Parameter_estimation/model_", model_name, "/")
  if(!dir.exists(dir_julia)){
    data_julia <- tibble()
    sprintf("Julia benchmark results for %s does not exist", model_name)
  }else{
    data_julia <- read_csv(str_c(dir_julia, "Estimation_statistics.csv"), col_types = cols()) |> 
      mutate(Alg = case_when(Alg == "IpoptAutoHess" ~ "Ipopt full H", 
                             Alg == "IpoptBlockAutoHess" ~ "Ipopt block H", 
                             Alg == "optimIPNewtonAutoHess" ~ "Optim IPNewton full H", 
                             Alg == "optimIPNewtonBlockAutoHess" ~ "Optim IPNewton block H", 
                             Alg == "OptimIPNewtonGN" ~ "Optim IPNewton FIM", 
                             Alg == "optimLBFGS" ~ "Optim LBFGS", 
                             Alg == "FidesAutoHess" ~ "Fides full H", 
                             Alg == "FidesAutoHessBlock" ~ "Fides block H", 
                             Alg == "FidesBFGS" ~ "Fides BFGS", 
                             Alg == "FidesGN" ~ "Fides FIM"))  |> 
      mutate(software = "Julia")
  }

  if(!is.null(jl_solver)) data_julia = data_julia |> filter(Solver == jl_solver)
  
  data_pyPesto <- get_pyPesto_results(model_name) |> 
    mutate(software = "PyPesto")
  
  data_tot <- get_rank(data_julia |> 
                         select(Cost, Run_time, Alg, Solver, N_iter, software) |> 
                         bind_rows(data_pyPesto) |> 
                         mutate(Cost_norm = Cost - min(Cost) + 1))
  
  data_run_time <- data_tot |> group_by(Alg) |> summarise(median_run_time = median(Run_time, na.rm = T))
  pos <- data_run_time$Alg[order(data_run_time$median_run_time)]
  
  p1 <- ggplot(data_tot, aes(Alg, Run_time, fill=software)) + 
    geom_violin() + 
    geom_jitter(width=0.1, size=0.5) + 
    geom_crossbar(stat="summary", fun = median, width=0.5, color=col_highlight[5], linewidth=1.0) +
    scale_y_log10() + 
    scale_x_discrete(limits = pos) + 
    scale_fill_manual(values = cbPalette, name = "Software") +
    labs(x = "", y = "Run time [s]", title = str_c("Run time ", model_name)) +
    my_theme + 
    theme(plot.background = element_rect(fill = "white"))
  
  p2 <- ggplot(data_tot, aes(rank, Cost_norm, color = Alg)) + 
    geom_line(linewidth=1.5) + 
    geom_point(size=1.5) + 
    scale_color_manual(values = cbPalette, name = "Optimizer") + 
    scale_y_log10(limits=ylim) + 
    scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 100)) +
    labs(x = "", y = "Normalised likelihood", title = str_c("Waterfall ", model_name)) + 
    theme_bw(base_size = 16) + 
    theme(legend.position = "bottom") +
    theme(plot.background = element_rect(fill = "white"))
  
  data_cp <- data_tot |> rename("Alg_cp" = "Alg")
  p3 <- ggplot(data_tot, aes(Run_time, Cost_norm)) +
    geom_point(data=data_cp, mapping=aes(Run_time, Cost_norm), color = col_highlight[1]) +
    geom_point(color = col_highlight[5]) + 
    facet_wrap(~Alg) + 
    scale_x_log10() + 
    scale_y_log10(limits=ylim) + 
    labs(x = "Run time [s]", y = "Normalised cost", title = str_c("Cost and run time  ", model_name)) + 
    theme_bw(base_size = 16) +
    theme(plot.background = element_rect(fill = "white"))
  
  # Generate convergence data output as html table 
  data_convergence <- data_tot |> 
    mutate(has_converged = Cost < min(Cost) + 0.1) |> 
    filter(has_converged == TRUE) |> 
    group_by(Alg) |> 
    summarise(n_convergence = n(), 
              average_time = mean(Run_time), 
              average_iterations = mean(N_iter))
  my_table <- data_convergence |> 
    gt() |> 
    tab_header(title = str_c(model_name), "Converged starts") |> 
    fmt(columns = starts_with("average_time"), fns = function(x) sprintf("%.2e", x)) |> 
    fmt(columns = starts_with("average_iterations"), fns = function(x) sprintf("%.1f", x)) |> 
    cols_label(Alg= "Algorithm", n_convergence="Number of converged starts", average_time = "Time per converged start [s]")
  
  data_not_convergence <- data_tot |> 
    mutate(has_converged = Cost < min(Cost) + 0.1) |> 
    filter(has_converged == FALSE) |> 
    group_by(Alg) |> 
    summarise(n_convergence = n(), 
              average_time = mean(Run_time), 
              average_iterations = mean(N_iter))
  my_table_not <- data_not_convergence |> 
    gt() |> 
    tab_header(title = str_c(model_name, " : None converged starts")) |> 
    fmt(columns = starts_with("average_time"), fns = function(x) sprintf("%.2e", x)) |> 
    fmt(columns = starts_with("average_iterations"), fns = function(x) sprintf("%.1f", x)) |> 
    cols_label(Alg= "Algorithm", n_convergence="Number of non converged starts", average_time = "Time per non converged start [s]")
  
  dir_save <- str_c("./Results/Parameter_estimation/", model_name, "/")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  
  ggsave(str_c(dir_save, "Run_time.png"), p1, width = BASE_WIDTH*2.3, height = BASE_HEIGHT*2, dpi=300)
  ggsave(str_c(dir_save, "Waterfall.png"), p2, width = BASE_WIDTH*2, height = BASE_HEIGHT*2, dpi=300)
  ggsave(str_c(dir_save, "Cost_vs_time.png"), p3, width = BASE_WIDTH*2, height = BASE_HEIGHT*2, dpi=300)
  gtsave(my_table, str_c(dir_save, "Convergence_statistics.html"))
  gtsave(my_table_not, str_c(dir_save, "Non_convergence_statistics.html"))
}


process_model("Boehm_JProteomeRes2014")

process_model("Fiedler_BMC2016", ylim=c(0.9, 1e3), jl_solver="QNDF")

process_model("Fujita_SciSignal2010", ylim=c(0.9, 1e9))

