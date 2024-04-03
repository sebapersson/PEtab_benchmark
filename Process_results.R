library(tidyverse)
library(ggthemes)
library(gt)
library(rhdf5)
library(cowplot)
library(ggrepel)


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

get_pyPesto_results <- function(model_name)
{
  
  pyPestoFiles <- list.files("pypesto_benchmark/results")[str_detect(list.files("pypesto_benchmark/results"), str_c(model_name, ".+", ".csv"))]
  pyPestoStatFiles <- str_c("pypesto_benchmark/stats/",  list.files("pypesto_benchmark/stats")[str_detect(list.files("pypesto_benchmark/stats"), str_c(model_name, ".+", ".hdf5"))])
  pyPestoData <- tibble()
  for(i in 1:length(pyPestoFiles)){
  
    pyPestoData_ <- read_csv(str_c("pypesto_benchmark/results/", pyPestoFiles[i]), col_names = F, col_types = cols()) |> 
      filter(!is.infinite(X1))
    if(length(names(pyPestoData_)) == 2){
      names(pyPestoData_) <- c("Cost", "Run_time")
    }else{
      names(pyPestoData_) <- c("Cost", "Run_time", "Cost_x0")      
    }
      
    pyPestoStatFile <- h5read(pyPestoStatFiles[i], "/")
    pyPestoStatsNames <- names(pyPestoStatFile)
    pyPestoStatsNames <- pyPestoStatsNames[get_name_order(pyPestoStatFile, pyPestoStatsNames)]
    pyPestoData_ <- pyPestoData_[1:length(pyPestoStatsNames), ]

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
                             Alg == "optimIPNewtonAutoHess" ~ "Optim full H", 
                             Alg == "optimIPNewtonBlockAutoHess" ~ "Optim block H", 
                             Alg == "OptimIPNewtonGN" ~ "Optim FIM", 
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
  
  data_tot <- get_rank(drop_na(data_julia) |> 
                         select(Cost, Run_time, Alg, Solver, N_iter, software) |> 
                         bind_rows(data_pyPesto) |> 
                         mutate(Cost_norm = Cost - min(Cost) + 1)) |> 
    filter(!is.infinite(Cost)) |> 
    filter(Run_time != 0) |> 
    mutate(model_name = model_name)
  
  data_run_time <- data_tot |> group_by(Alg, software) |> summarise(median_run_time = median(Run_time, na.rm = T))
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
  
  # Normalized run-time
  data_run_time <- data_run_time |> 
    mutate(n_time = median_run_time / min(data_run_time$median_run_time)) |> 
    mutate(str_write = sprintf("%.1f", n_time))
    
  p1_ <- ggplot(data_run_time, aes(Alg, n_time, fill = software)) + 
    geom_bar(stat="identity") + 
    geom_text(aes(y = n_time+0.3, label=str_write), size=8) +
    scale_x_discrete(limits = pos) + 
    labs(x = "", y = "Median run-time normalized by fastest algorithm", title = str_c("Normalised median run time ", model_name)) +
    scale_fill_manual(values = cbPalette[-1], name = "Software") + 
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
  ggsave(str_c(dir_save, "Run_time_norm.png"), p1_, width = BASE_WIDTH*2.3, height = BASE_HEIGHT*2, dpi=300)
  ggsave(str_c(dir_save, "Waterfall.png"), p2, width = BASE_WIDTH*1.5, height = BASE_HEIGHT*2, dpi=300)
  ggsave(str_c(dir_save, "Cost_vs_time.png"), p3, width = BASE_WIDTH*2, height = BASE_HEIGHT*2, dpi=300)
  gtsave(my_table, str_c(dir_save, "Convergence_statistics.html"))
  gtsave(my_table, str_c(dir_save, "Convergence_statistics.png"))
  gtsave(my_table_not, str_c(dir_save, "Non_convergence_statistics.html"))
  write_csv(data_tot, str_c(dir_save, "Processed_data.csv"))
  
  # SVG for paper 
  ggsave(str_c(dir_save, "Run_time.svg"), p1, width = BASE_WIDTH*2.3, height = BASE_HEIGHT*2, dpi=300)
  ggsave(str_c(dir_save, "Run_time_norm.svg"), p1_, width = BASE_WIDTH*2.3, height = BASE_HEIGHT*2, dpi=300)
  ggsave(str_c(dir_save, "Waterfall.svg"), p2, width = BASE_WIDTH*2, height = BASE_WIDTH*2, dpi=300)
  ggsave(str_c(dir_save, "Cost_vs_time.svg"), p3, width = BASE_WIDTH*2, height = BASE_HEIGHT*2, dpi=300)
  
  p_table <- ggdraw() + draw_image(str_c(dir_save, "Convergence_statistics.png"), scale = 0.8)
  p_summary <- ggpubr::ggarrange(p1, p2, p3, p_table, widths = c(1.3, 1))
  ggsave(str_c(dir_save, "Summary", model_name, ".png"), p_summary, width = BASE_WIDTH*3.5, height = BASE_HEIGHT*3.5, dpi=300, bg="white")
  
}


process_model("Boehm_JProteomeRes2014")

process_model("Fiedler_BMC2016", ylim=c(0.9, 1e3), jl_solver="QNDF")

process_model("Fujita_SciSignal2010", ylim=c(0.9, 1e9))

process_model("Brannmark_JBC2010", ylim=c(0.9, 1e9))

process_model("Bruno_JExpBot2016")

process_model("Weber_BMC2015", ylim=c(0.9, 1e9))

process_model("Zheng_PNAS2012", ylim=c(0.9, 1e5))

process_model("Beer_MolBioSystems2014")

process_model("Sneyd_PNAS2002")

process_model("Bachmann_MSB2011")

process_model("Crauste_CellSystems2017")

process_model("Schwen_PONE2014")

process_model("Lucarelli_CellSystems2018")

process_model("Borghans_BiophysChem1997")

process_model("Okuonghae_ChaosSolitonsFractals2020")

process_model("Rahman_MBS2016")

process_model("Isensee_JCB2018")

model_name <- "Isensee_JCB2018"

process_model("Oliveira_NatCommun2021")


# ---------------------------------------------------------------------------------------------------------
# Process subproblem solver for fides 
# ---------------------------------------------------------------------------------------------------------
model_names = c("Boehm_JProteomeRes2014", "Fiedler_BMC2016", "Crauste_CellSystems2017", "Fujita_SciSignal2010")
data_plot = tibble()
for(model_name in model_names){
  dir_julia <- str_c("Master-Thesis/Intermediate/Benchmarks/Parameter_estimation/model_", model_name, "/")
  if(!dir.exists(dir_julia)){
    data_julia <- tibble()
    sprintf("Julia benchmark results for %s does not exist", model_name)
  }else{
    data_julia <- read_csv(str_c(dir_julia, "Estimation_statisticsnot_full.csv"), col_types = cols()) |> 
      mutate(Alg = case_when(Alg == "FidesAutoHess" ~ "Fides full H", 
                             Alg == "FidesGN" ~ "Fides FIM"), 
             model_name = model_name)  |> 
      mutate(software = "Julia", 
             fullsubproblem = F) |> 
      select(-Ret_code)
    data_julia_full = read_csv(str_c(dir_julia, "Estimation_statisticsfull.csv"), col_types = cols()) |> 
      mutate(Alg = case_when(Alg == "FidesAutoHess" ~ "Fides full H", 
                             Alg == "FidesGN" ~ "Fides FIM"), 
             model_name = model_name)  |> 
      mutate(software = "Julia", 
             fullsubproblem = T) |> 
      select(-Ret_code)
  }
  data_plot = bind_rows(data_plot, data_julia) |> 
    bind_rows(data_julia_full) |> 
    filter(Alg %in% c("Fides full H", "Fides FIM")) |> 
    filter(!is.infinite(Run_time))  
}


min_cost = min(data_plot$Cost, na.rm = T)
data_converged = data_plot |> 
  group_by(Alg, fullsubproblem, model_name) |> 
  summarise(n_converged = sum(min(Cost, na.rm = T) + 0.1 > Cost), 
            runtime = mean(Run_time, na.rm=T)*1000) 

p1 = ggplot(data_converged, aes(model_name, n_converged, fill = fullsubproblem)) + 
  geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~Alg, ncol = 2) + 
  scale_fill_manual(values = cbPalette[-1], name = "Full subproblem") +
  scale_y_continuous(expand=c(0, 0)) + 
  labs(x = "", y = "Number of converged starts", title = "Converged starts (1000 multistarts)") + 
  theme_bw(base_size = 18)

p2 = ggplot(data_plot, aes(model_name, Run_time, fill = fullsubproblem)) + 
  geom_violin(draw_quantiles = 0.5) + 
  facet_wrap(~Alg, ncol = 2) + 
  scale_fill_manual(values = cbPalette[-1], name = "Full subproblem") +
  scale_y_log10() + 
  labs(x = "", y = "Run time [s]", title = "Runtime (1000 multistarts)") + 
  theme_bw(base_size = 18)

dir_save <- "Results/Parameter_estimation/Summary/"
if(!dir.exists(dir_save)) dir.create(dir_save)
ggsave(str_c(dir_save, "Fides_subproblem.png"), p1, dpi=300, height = BASE_HEIGHT*2, width = BASE_WIDTH*3.5)
ggsave(str_c(dir_save, "Fides_subproblem_time.png"), p2, dpi=300, height = BASE_HEIGHT*2, width = BASE_WIDTH*3.5)

# ---------------------------------------------------------------------------------------------------------
# Process across all models 
# ---------------------------------------------------------------------------------------------------------
dir_models <- str_c("Results/Parameter_estimation/", list.files("Results/Parameter_estimation/")) 
data_tot <- tibble()
for(i in 1:length(dir_models)){
  if(dir_models[i] == "Results/Parameter_estimation/Summary") next
  file_read <- str_c(dir_models[i], "/Processed_data.csv")
  data_tot <- bind_rows(data_tot, read_csv(file_read, col_types = cols()))
}
data_tot <- filter(data_tot, Alg != "Optim block H")

# Compute convergence statistics per model 
model_list <- unique(data_tot$model_name)
data_converged <- tibble()
for(i in 1:length(model_list)){
  data_tmp <- data_tot |> 
    filter(model_name == model_list[i]) 
  min_value = min(data_tmp$Cost, na.rm = T)
  data_tmp <- data_tmp |> 
    group_by(Alg, Solver, software, model_name) |> 
    summarise(n_converged_starts = sum(Cost < min_value + 0.1), 
              avg_iter = sum((Cost < min_value + 0.1) * N_iter / sum(Cost < min_value + 0.1)), 
              avg_time = mean(Run_time, na.rm=T)) |> 
    mutate(t_succ = avg_time / n_converged_starts) #|>
  min_t_succ <- min(data_tmp$t_succ, na.rm = T)
  max_conv = max(data_tmp$n_converged_starts, na.rm = T)
  max_time <- max(data_tmp$avg_time, na.rm = T)
  min_time <- min(data_tmp$avg_time, na.rm = T)
  data_tmp <- mutate(data_tmp, 
                     OE = min_t_succ / t_succ, 
                     conv_ratio = n_converged_starts / max_conv, 
                     time_ratio = avg_time / min_time)
  data_converged <- bind_rows(data_converged, data_tmp)
}

# --------------------------------------
# Aggregate converged results
# --------------------------------------
pos <- c("Bachmann_MSB2011", "Beer_MolBioSystems2014", "Boehm_JProteomeRes2014", "Brannmark_JBC2010", "Bruno_JExpBot2016", "Fiedler_BMC2016", 
         "Fujita_SciSignal2010", "Sneyd_PNAS2002", "Weber_BMC2015", "Zheng_PNAS2012", "Crauste_CellSystems2017", "Elowitz_Nature2000", 
         "Isensee_JCB2018", "Rahman_MBS2016", "Borghans_BiophysChem1997", "Okuonghae_ChaosSolFract2020", "Schwen_PONE2014", 
         "Lucarelli_CellSystems2018", "Oliveira_NatCommun2021")
pos_short <- str_sub(str_extract(pos, "[:alpha:]+_"), 1, -2)
n_param <- c(113, 72, 9, 22, 13, 22, 19, 15, 36, 46, 12, 21, 46, 7, 23, 16, 30, 84, 12)
pos <- pos[order(n_param)]
pos_short <- pos_short[order(n_param)]

data_converged <- data_converged |> 
  complete(model_name, Alg) |> 
  mutate(software = case_when(Alg %in% c("PyPesto Fides BFGS", "PyPesto Fides FIM", "PyPesto Fides") ~ "PyPesto", 
                              T ~ "Julia"))
data_converged$model_short <- str_sub(str_extract(data_converged$model_name, "[:alpha:]+_"), 1, -2)

p1 <- ggplot(data_converged, aes(model_short, Alg, fill = OE)) + 
  geom_tile(color="white") + 
  geom_text(aes(label = ifelse(is.na(OE), "", sprintf("%.2f", OE))), size=5.0, color="#CC333F") +
  scale_fill_viridis_c() + 
  facet_grid(rows=vars(software), scales = "free_y", space="free_y") +
  scale_x_discrete(limits = pos_short, expand=c(0, 0)) +  
  scale_y_discrete(expand=c(0, 0)) +
  labs(x = "", y = "") +
  theme_bw(base_size = 20.0) + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle=90, size=25), 
        axis.text.y = element_text(size=25))

p2 <- ggplot(data_converged, aes(model_short, Alg, fill = conv_ratio)) + 
  geom_tile(color="white") + 
  geom_text(aes(label = ifelse(is.na(OE), "", sprintf("%d", n_converged_starts))), size=5.0, color="#CC333F") +
  scale_fill_viridis_c(option="C") + 
  facet_grid(rows=vars(software), scales = "free_y", space="free_y") +
  scale_x_discrete(limits = pos_short, expand=c(0, 0)) +  
  scale_y_discrete(expand=c(0, 0)) +
  labs(x = "", y = "") +
  theme_bw(base_size = 20.0) + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle=90, size=25), 
        axis.text.y = element_text(size=25))

p3 <- ggplot(data_converged, aes(model_short, Alg, fill = time_ratio)) + 
  geom_tile(color="white") + 
  geom_text(aes(label = ifelse(is.na(OE), "", sprintf("%.2f", avg_time*1000 / 3600))), size=5.0, color="#CC2A41") +
  scale_fill_viridis_c(option="E", direction=-1, trans="log10") + 
  facet_grid(rows=vars(software), scales = "free_y", space="free_y") +
  scale_x_discrete(limits = pos_short, expand=c(0, 0)) +  
  scale_y_discrete(expand=c(0, 0)) +
  labs(x = "", y = "") +
  theme_bw(base_size = 20.0) + 
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle=90, size=25), 
        axis.text.y = element_text(size=25))

ggsave(str_c(dir_save, "OE_index.svg"), p1, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*2.5)
ggsave(str_c(dir_save, "Converged.svg"), p2, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*2.5)
ggsave(str_c(dir_save, "Run_time.svg"), p3, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*2.5)

# --------------------------------------
# Direct comparison using fides 
# --------------------------------------
data_bfgs_fim <- data_tot |> 
  filter(Alg %in% c("Fides BFGS", "Fides FIM", "PyPesto Fides BFGS", "PyPesto Fides FIM")) |> 
  group_by(model_name, Alg, software) |> 
  summarise(median_run_time = median(Run_time, na.rm = T), 
            mean_run_time = mean(Run_time, na.rm = T),
            mad_run_time = mad(Run_time, na.rm = T))

pyPesto_data <- data_bfgs_fim |> 
  filter(Alg %in% c("PyPesto Fides BFGS", "PyPesto Fides FIM")) |> 
  mutate(Alg = case_when(Alg == "PyPesto Fides BFGS" ~ "Fides BFGS", 
                         Alg == "PyPesto Fides FIM" ~ "Fides FIM")) |> 
  rename("median_run_time_pyPesto" = "median_run_time", 
         "mean_run_time_pyPesto" = "mean_run_time",
         "mad_run_time_pyPesto" = "mad_run_time") |> 
  select(-software)
julia_data <- data_bfgs_fim |> 
  filter(Alg %in% c("Fides BFGS", "Fides FIM"))

data_plot <- inner_join(julia_data, pyPesto_data, by = c("model_name", "Alg"))
max_val <- 5e3
data_line <- tibble(x = seq(from = 0.1, to = max_val, by = 1.0), y = (x*2.0), label="1") |> 
  bind_rows(tibble(x = seq(from = 0.1, to = max_val, by = 1.0), y = (x*4.0), label="2")) |> 
  bind_rows(tibble(x = seq(from = 0.01, to = max_val, by = 0.1), y = (x*0.5), label="3")) 

range_data = tibble(median_run_time=c(.1, max_val), median_run_time_pyPesto=c(.1, max_val))
p1 <- ggplot(data_plot, aes(x=mean_run_time, y=mean_run_time_pyPesto)) + 
  geom_abline(intercept = 0, slope = 1.0, linewidth=1.5, color=cbPalette[1]) + 
  geom_line(data=data_line, mapping=aes(x=x, y=y, group=label), linetype=2, linewidth=1.25, color=cbPalette[1]) +
  geom_point(aes(color = Alg), size=7.0) + 
  geom_text_repel(aes(label = model_name),  max.overlaps=1000, nudge_x=0.15) + 
  scale_x_log10(limits=c(0.1, max_val), expand=c(0, 0, 0.01, 0.1), breaks=c(1, 10, 100, max_val)) + 
  scale_y_log10(limits=c(0.1, max_val), expand=c(0, 0, 0.01, 0.1), breaks=c(1, 10, 100, max_val)) + 
  labs(x = "Average run time PEtab.jl [s]", y = "Average run time pyPesto [s]") +
  scale_color_manual(values = cbPalette[-1]) + 
  geom_rangeframe(data=range_data, mapping=aes(median_run_time, median_run_time_pyPesto), size=1.5) +
  my_theme + theme(panel.grid.major.y = element_blank(), legend.position = "none")
data_sum <- data_plot |> 
  group_by(Alg) |> 
  summarise(n_win = sum(mean_run_time < mean_run_time_pyPesto), 
            n_models = n(),
            speedup = mean(mean_run_time_pyPesto / mean_run_time))

# Convergence data 
data_bfgs_fim_conv <- data_converged |> 
  filter(Alg %in% c("Fides BFGS", "Fides FIM", "PyPesto Fides BFGS", "PyPesto Fides FIM")) |> 
  group_by(model_name, Alg, software) 

pyPesto_data_conv <- data_bfgs_fim_conv |> 
  filter(Alg %in% c("PyPesto Fides BFGS", "PyPesto Fides FIM")) |> 
  mutate(Alg = case_when(Alg == "PyPesto Fides BFGS" ~ "Fides BFGS", 
                         Alg == "PyPesto Fides FIM" ~ "Fides FIM")) |> 
  rename("avg_time_pyPesto" = "avg_time", 
         "n_converged_starts_pyPesto" = "n_converged_starts", 
         "avg_iter_pyPesto" = "avg_iter") 
julia_data_conv <- data_bfgs_fim_conv |> 
  filter(Alg %in% c("Fides BFGS", "Fides FIM"))

data_plot_conv <- inner_join(julia_data_conv, pyPesto_data_conv, by = c("model_name", "Alg"))
max_val <- 1100
range_data = tibble(x=0:max_val, y=0:max_val)
p2 <- ggplot(data_plot_conv, aes(x=n_converged_starts, y=n_converged_starts_pyPesto)) + 
  geom_abline(intercept = 0, slope = 1.0, linewidth=1.5, color=cbPalette[1]) + 
  geom_point(aes(color = Alg), size=7.0) + 
  scale_color_manual(values = cbPalette[-1]) + 
  scale_x_sqrt(limits = c(0, max_val), breaks=seq(from = 0, by = 100, to = max_val)) + 
  scale_y_sqrt(limits = c(0, max_val), breaks=seq(from = 0, by = 100, to = max_val)) + 
  labs(x = "Number of converged starts PEtab.jl", y = "Number of converged starts pyPesto") +
  geom_rangeframe(data=range_data, mapping=aes(x, y), size=1.5) +
  my_theme + theme(panel.grid.major.y = element_blank(), legend.position = "none")

max_val <- 3000
max(data_plot_conv$n_converged_starts_pyPesto / data_plot_conv$avg_time_pyPesto, na.rm = T)
p3 = ggplot(data_plot_conv, aes(y=n_converged_starts / avg_time + 0.1, x=n_converged_starts_pyPesto / avg_time_pyPesto + 0.1)) + 
  geom_abline(intercept = 0, slope = 1.0, linewidth=1.5, color=cbPalette[1]) + 
  geom_line(data=data_line, mapping=aes(x=x, y=y, group=label), linetype=2, linewidth=1.25, color=cbPalette[1]) +
  geom_point(aes(color = Alg), size=7.0) + 
  geom_text_repel(aes(label = model_name),  max.overlaps=1000, nudge_x=0.15) + 
  scale_x_log10(limits=c(0.05, max_val), expand=c(0, 0, 0.01, 0.1), breaks=c(1, 10, max_val, 1)) + 
  scale_y_log10(limits=c(0.05, max_val), expand=c(0, 0, 0.01, 0.1), breaks=c(1, 10, max_val, 1)) + 
  labs(x = "Converged starts / avg. run-time PEtab.jl", y = "Converged starts / avg. run-time pyPesto") +
  scale_color_manual(values = cbPalette[-1]) + 
  my_theme + theme(panel.grid.major.y = element_blank(), legend.position = "none")

ggsave(str_c(dir_save, "Fides_comparison_runtime.svg"), p1, width = BASE_WIDTH*1.5, height = BASE_WIDTH*1.5)
ggsave(str_c(dir_save, "Fides_comparison_cost.svg"), p2, width = BASE_WIDTH*1.5, height = BASE_WIDTH*1.5)
ggsave(str_c(dir_save, "Fides_comparison_OE.svg"), p3, width = BASE_WIDTH*1.5, height = BASE_WIDTH*1.5)

# --------------------------------------------------------------
# Direct pairwise comparison of different optimizer 
# --------------------------------------------------------------
compare_opt_pair <- function(opt1, opt2, data_tot, add_time=1000, log_axis=F)
{
  
  data_opt1 <- data_converged |> 
    filter(Alg == opt1) |> 
    select(-Solver, -software, -OE, -conv_ratio, -time_ratio, -model_short) |>
    rename("n_converged_starts_opt1" = "n_converged_starts", 
           "avg_iter_opt1" = "avg_iter", 
           "t_succ_opt1" = "t_succ", 
           "avg_time_opt1" = "avg_time")
  data_opt2 <- data_converged |> 
    filter(Alg == opt2) |> 
    select(-Solver, -software, -OE, -conv_ratio, -time_ratio, -model_short) |>
    rename("n_converged_starts_opt2" = "n_converged_starts", 
           "avg_iter_opt2" = "avg_iter", 
           "t_succ_opt2" = "t_succ", 
           "avg_time_opt2" = "avg_time")
  
  data_plot <- inner_join(data_opt1, data_opt2, by = c("model_name"))
  
  max_val <- max(max(data_plot$n_converged_starts_opt1, na.rm = T), max(data_plot$n_converged_starts_opt2, na.rm = T))
  range_data <- tibble(x = c(0, max_val), y = c(0, max_val))
  p1 <- ggplot(data_plot, aes(n_converged_starts_opt1, n_converged_starts_opt2)) + 
    geom_point(size=3.0) + 
    geom_abline(intercept = 0, slope = 1.0, linewidth=1.5, color=cbPalette[1]) + 
    labs(x = opt1, y = opt2, 
         title = str_c(opt1, " vs ", opt2, " : Converged starts")) + 
    geom_rangeframe(data=range_data, mapping=aes(x, y), size=1.5) +
    geom_text_repel(aes(label = model_name),  max.overlaps=1000, nudge_x=0.25) + 
    my_theme + theme(panel.grid.major.y = element_blank(), legend.position = "none")
  if(log_axis == T){
    p1 <- p1 + 
      scale_x_sqrt(limits=c(0.0, max_val)) + 
      scale_y_sqrt(limits=c(0.0, max_val)) 
  }
  
  max_val <- max(max(data_plot$avg_time_opt1, na.rm = T), max(data_plot$avg_time_opt2, na.rm = T))
  data_line <- tibble(x = seq(from = 0.01, to = max_val+add_time, by = 0.1), y = (x*2.0), label="1") |> 
    bind_rows(tibble(x = seq(from = 0.01, to = max_val+add_time, by = 0.1), y = (x*0.5), label="2"))
  range_data <- tibble(x = c(0, max_val+add_time), y = c(0, max_val+add_time))
  p2 <- ggplot(data_plot, aes(avg_time_opt1, avg_time_opt2)) + 
    geom_point(size=4.0) + 
    geom_abline(intercept = 0, slope = 1.0, linewidth=1.5, color=cbPalette[1]) + 
    labs(x = opt1, y = opt2) + 
    scale_x_log10(limits=c(1.0, max_val+add_time)) + 
    scale_y_log10(limits=c(1.0, max_val+add_time)) + 
    geom_line(data=data_line, mapping=aes(x=x, y=y, group=label), linetype=2, linewidth=1.25, color=cbPalette[1]) +
    geom_rangeframe(data=range_data, mapping=aes(x, y), size=1.5) +
    geom_text_repel(aes(label = model_name),  max.overlaps=1000, nudge_x=0.15) + 
    labs(x = str_c(opt1, " avg. run time [s]"), y = str_c(opt2, " avg. run time [s]"), 
         title = str_c(opt1, " vs ", opt2, " : Run time")) +
    my_theme + theme(panel.grid.major.y = element_blank(), legend.position = "none")
  
  dir_save <- "Results/Parameter_estimation/Summary/Compare_opt_pair/"
  if(!dir.exists(dir_save)) dir.create(dir_save)
  ggsave(str_c(dir_save, "Compare_conv_", opt1, "_", opt2, ".svg"), p1, width = BASE_WIDTH*1.5, height = BASE_WIDTH*1.5)
  ggsave(str_c(dir_save, "Compare_time_", opt1, "_", opt2, ".svg"), p2, width = BASE_WIDTH*1.5, height = BASE_WIDTH*1.5)
  return(data_plot)
}


data1 <- compare_opt_pair("Fides full H", "Optim full H", data_tot, add_time=1000, log_axis=T)
data2 <- compare_opt_pair("Fides FIM", "Optim FIM", data_tot, add_time=1000, log_axis=T)
data3 <- compare_opt_pair("PyPesto Fides FIM", "Optim FIM", data_tot, add_time=1000, log_axis=T)

data1$full_h <- "Yes"
data3$full_h <- "No"
data_plot <- bind_rows(data1, data3)
max_val <- 1100
data_line <- tibble(x = c(0, max_val), 
                    y = c(0, max_val))

p <- ggplot(data_plot, aes(n_converged_starts_opt1, n_converged_starts_opt2)) + 
  geom_point(aes(color=full_h), size=7.0) + 
  geom_line(data=data_line, mapping=aes(x=x, y=y), color=cbPalette[1], linewidth=2.5) + 
  #geom_text_repel(aes(label = model_name),  max.overlaps=1000, nudge_x=0.25, color="black") + 
  scale_x_sqrt(limits=c(0, max_val), name = "Number of converged starts Newton trust-region", breaks=seq(from=0, by = 100, to=1000)) + 
  scale_y_sqrt(limits=c(0, max_val), name = "Number of converged starts Newton interior point", breaks=seq(from=0, by = 100, to=1000)) +
  scale_color_manual(values = col_highlight[c(3, 5)]) +
  theme_classic(base_size = 20) + 
  theme(legend.position = "none", 
        panel.grid.major=element_blank()) 

data_sum <- data_plot |> 
  group_by(Alg.x) |> 
  summarise(n_best_int = sum(n_converged_starts_opt1 > n_converged_starts_opt2, na.rm = T), 
            n_equal = sum(n_converged_starts_opt1 == n_converged_starts_opt2, na.rm = T),
            n = n())

data_fim <- data_plot |> 
  filter(Alg.x != "Fides full H") |> 

ggsave(str_c(dir_save, "Trust_vs_interior.svg"), p, width = BASE_WIDTH*1.5, height = BASE_WIDTH*1.5)
