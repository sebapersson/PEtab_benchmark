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


# ---------------------------------------------------------------------------------------------------------------
# Run-times Julia and AMICI 
# ---------------------------------------------------------------------------------------------------------------
to_wide_data = function(data)
{
  data_amici = data |> filter(solver == "AMICI") |> 
    rename("Time_amici" = "Time") |> 
    select(Model, Time_amici, n_parameters)
  data_QNDF = data |> filter(solver == "QNDF") |> 
    select(Model, Time, n_parameters)
  data_Rodas = data |> filter(solver == "Rodas5P") |> 
    select(Model, Time, n_parameters)
  data_p1 = inner_join(data_QNDF, data_amici, by = c("Model", "n_parameters")) |> 
    mutate(solver = "QNDF")
  data_p2 = inner_join(data_Rodas, data_amici, by = c("Model", "n_parameters")) |> 
    mutate(solver = "Rodas5P")
  return(bind_rows(data_p1, data_p2))
}

# Read Julia data 
dir_data <- "./Master-Thesis/Intermediate/Benchmarks/Cost_grad_hess/"

data_cost_gradient <- read_csv(str_c(dir_data, "Gradient_cost_small_models.csv"), col_types = cols()) |> 
  group_by(Method_info, Model, absTol, relTol, solver) |> 
  summarise(Time = median(Time)) |> 
  mutate(Model = str_sub(Model, start=7)) |> 
  filter(solver != "Vern7(Rodas5P)") |> 
  filter(solver != "KenCarp47") |> 
  filter(Model != "Giordano_Nature2020")
  
data_cost <- data_cost_gradient |> filter(Method_info == "Standard") |> ungroup() |>  select(-Method_info)
data_gradient <- data_cost_gradient |> filter(Method_info == "ForwardDiff") |> ungroup() |> select(-Method_info)
data_gradient_forward <- data_cost_gradient |> filter(Method_info == "ForEq_AutoDiff") |> ungroup() |> select(-Method_info)


# Read PyPesto data 
dir_pypesto <- "./Intermediate/AMICI_cost_grad/"
path_files <- str_c(dir_pypesto, list.files(dir_pypesto))
data_pypesto <- tibble()
for(file in path_files) data_pypesto <- bind_rows(data_pypesto, read_csv(file, col_types = cols()))
names(data_pypesto) <- c("Model", "Cost", "Forward", "Adj", "n_parameters")
data_pypesto <- data_pypesto |> 
  mutate(absTol=1e-8, relTol=1e-8, solver="AMICI") |> 
  group_by(Model, absTol, relTol, solver) |> 
  summarise(Cost = median(Cost), 
         Forward = median(Forward), 
         Adj = median(Adj), 
         n_parameters = median(n_parameters)) |> 
  filter(Model != "Giordano_Nature2020") |> 
  ungroup()

data_pypesto_cost <- data_pypesto |> select(-c("Forward", "Adj")) |> rename("Cost_amici" = "Cost")
data_pypesto_grad <- data_pypesto |> select(-c("Cost", "Adj")) |> rename("Gradient_amici" = "Forward")
data_pypesto_grad_adj <- data_pypesto |> select(-c("Forward", "Cost")) |> rename("Gradient_amici_adj" = "Adj")

data_pypesto_cost_ <- data_pypesto_cost |> select(-n_parameters) |> rename("Time" = "Cost_amici")
data_pypesto_grad_ <- data_pypesto_grad |> select(-n_parameters) |> rename("Time" = "Gradient_amici")
data_pypesto_grad_adj_ <- data_pypesto_grad_adj |> select(-n_parameters) |> rename("Time" = "Gradient_amici_adj")

data_plot_cost <- bind_rows(data_cost, data_pypesto_cost_) |> 
  inner_join(select(data_pypesto_cost, -solver), by=c("Model", "absTol", "relTol")) |> 
  filter(!is.infinite(Time))

data_plot_gradient <- bind_rows(data_gradient, data_pypesto_grad_) |> 
  inner_join(select(data_pypesto_grad, -solver), by=c("Model", "absTol", "relTol")) |> 
  filter(!is.infinite(Time))

data_plot_gradient_forward <- bind_rows(data_gradient_forward, data_pypesto_grad_) |> 
  inner_join(select(data_pypesto_grad, -solver), by=c("Model", "absTol", "relTol")) |> 
  filter(!is.infinite(Time))

data_plot_cost_alt = to_wide_data(data_plot_cost)
data_plot_gradient_forward_alt = to_wide_data(data_plot_gradient_forward)

pos <- unique(data_plot_cost$Model[order(data_plot_cost$n_parameters)])

p1 <- ggplot(data_plot_cost, aes(Model, Time / Cost_amici, fill = solver)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = cbPalette[-1], name = "ODE Solver") +
  scale_x_discrete(limits=pos, name="") + 
  scale_y_continuous(expand=c(0, 0, 0.01, 0.1), name = "", breaks=1:5) +
  geom_hline(yintercept = 1.0) + 
  coord_flip() + 
  my_minimal + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        legend.position = "none")

max_val = max(c(data_plot_cost_alt$Time, data_plot_cost_alt$Time_amici)) * 1.2
min_val = 1e-4
line1 = tibble(x=c(min_val, max_val), y=c(min_val, max_val))
line2 = tibble(x=seq(from=min_val, to=max_val, length.out=40), y=x*2.0)
line3 = tibble(x=seq(from=min_val, to=max_val, length.out=80), y=x*0.5)
p1_alt = ggplot(data_plot_cost_alt) + 
  geom_line(data=line1, mapping=aes(x=x, y=y), linewidth=1.0, color=cbPalette[1]) +
  geom_line(data=line2, mapping=aes(x=x, y=y), linewidth=1.0, color=cbPalette[1], linetype="dashed") +
  geom_line(data=line3, mapping=aes(x=x, y=y), linewidth=1.0, color=cbPalette[1], linetype="dashed") + 
  geom_point(aes(Time, Time_amici, color=solver), size = 4.0) + 
  scale_x_log10(limits = c(min_val, max_val), expand=c(0, 0)) + 
  scale_y_log10(limits = c(min_val, max_val), expand=c(0, 0)) + 
  scale_color_manual(values = cbPalette[-1]) + 
  labs(x = "Run time AMICI [s]", y = "Run time PEtab.jl [s]") + 
  theme_classic(base_size = 18)

data_sum_cost <- data_plot_cost_alt |> 
  mutate(is_small = n_parameters <= 20) |> 
  group_by(solver, is_small) |> 
  summarise(mean_time = median(Time), 
            mean_AMICI = median(Time_amici), 
            frac_best = sum(Time < Time_amici) / n(), 
            n = sum(is_small)) |> 
  mutate(ratio = mean_time / mean_AMICI)


p2 <- ggplot(data_plot_gradient, aes(Model, Time / Gradient_amici, fill = solver)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = cbPalette[-1], name = "ODE Solver") +
  scale_x_discrete(limits=pos, name="") + 
  scale_y_continuous(expand=c(0, 0, 0.01, 0.1), name = "", breaks=1:4) +
  geom_hline(yintercept = 1.0) + 
  coord_flip() + 
  my_minimal +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = "none")  

max_val = max(c(data_plot_gradient_forward_alt$Time, data_plot_gradient_forward_alt$Time_amici)) * 1.2
min_val = 1e-3
line1 = tibble(x=c(min_val, max_val), y=c(min_val, max_val))
line2 = tibble(x=seq(from=min_val, to=max_val, length.out=40), y=x*2.0)
line3 = tibble(x=seq(from=min_val, to=max_val, length.out=80), y=x*0.5)
line4 = tibble(x=seq(from=min_val, to=max_val, length.out=80), y=x*3.0)
p2_alt = ggplot(data_plot_gradient_forward_alt) + 
  geom_line(data=line1, mapping=aes(x=x, y=y), linewidth=1.0, color=cbPalette[1]) +
  geom_line(data=line2, mapping=aes(x=x, y=y), linewidth=1.0, color=cbPalette[1], linetype="dashed") +
  geom_line(data=line3, mapping=aes(x=x, y=y), linewidth=1.0, color=cbPalette[1], linetype="dashed") + 
  geom_line(data=line4, mapping=aes(x=x, y=y), linewidth=1.0, color=cbPalette[1], linetype="dashed") + 
  geom_point(aes(Time, Time_amici, color=solver), size = 4.0) + 
  scale_x_log10(limits = c(min_val, max_val), expand=c(0, 0)) + 
  scale_y_log10(limits = c(min_val, max_val), expand=c(0, 0)) + 
  scale_color_manual(values = cbPalette[-1]) + 
  labs(y = "Run time AMICI [s]", x = "Run time PEtab.jl [s]") + 
  theme_classic(base_size = 18)
data_sum_grad <- data_plot_gradient_forward_alt |> 
  group_by(solver) |> 
  summarise(mean_time = median(Time), 
            mean_AMICI = median(Time_amici), 
            frac_best = sum(Time < Time_amici) / n(), 
            n_best = sum(Time < Time_amici)) |> 
  mutate(ratio = 1 / (mean_time / mean_AMICI))

data_plot_gradient_forward <- data_plot_gradient_forward |> filter(Time / Gradient_amici < 5)
p2_ <- ggplot(data_plot_gradient_forward, aes(Model, Time / Gradient_amici, fill = solver)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = cbPalette[-1], name = "ODE Solver") +
  scale_x_discrete(limits=pos, name="") + 
  scale_y_continuous(expand=c(0, 0, 0.01, 0.1), name = "", breaks=1:4) +
  geom_hline(yintercept = 1.0) + 
  coord_flip() + 
  my_minimal +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = "none")  

# Get corresponding plot but scale with run-times 
data_cost__ <- data_plot_cost |> 
  select(Model, solver, Time) |> 
  rename("Time_cost" = "Time")
data_plot_scaling <- data_plot_gradient |> 
  inner_join(data_cost__, by=c("Model", "solver"))

p3 <- ggplot(data_plot_scaling, aes(Model, Time / Time_cost, fill = solver)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = cbPalette[-1], name = "ODE Solver") +
  scale_x_discrete(limits=pos, name="") + 
  scale_y_continuous(expand=c(0, 0, 0.01, 0.1), name = "") +
  geom_text(aes(y = n_parameters, label=sprintf("%d", n_parameters))) + 
  coord_flip() + 
  my_minimal + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = "none")  

max_val = max(data_plot_scaling$n_parameters) + 3
min_val = -0.6
line1 = tibble(x = c(0, max_val), y = x)
p3_alt = ggplot(data_plot_scaling, aes(n_parameters, Time / Time_cost, color = solver)) + 
  geom_point(size=4.0) + 
  geom_line(data=line1, mapping=aes(x, y), color = cbPalette[1], linewidth=1.0) + 
  labs(x = "Number of parameters to take gradient on", y = "Run time gradient / Run time cost") + 
  scale_x_continuous(expand=c(0, 0), limits=c(min_val, max_val)) + 
  scale_y_continuous(expand=c(0, 0), limits=c(min_val, max_val)) + 
  geom_text_repel(aes(label = Model),  max.overlaps=1000, nudge_x=0.15) + 
  scale_color_manual(values = cbPalette[c(4, 2, 3)]) + 
  theme_classic(base_size = 18) + 
  theme(panel.grid.major = element_line(linewidth=0.05, color="grey10"))
scaling_sum <- data_plot_scaling |> 
  group_by(solver) |> 
  summarise(n_less = sum(Time / Time_cost < n_parameters), 
            ratio = n_less / 18 * 100)

data_cost_hessian <- read_csv(str_c(dir_data, "Hessian_cost_small_models.csv"), col_types = cols()) |> 
  group_by(Method_info, Model, absTol, relTol, solver, What_calc) |> 
  summarise(Time = median(Time)) |> 
  ungroup() |> 
  mutate(Model = str_sub(Model, start=7)) |> 
  filter(solver != "Vern7(Rodas5P)") |> 
  filter(solver != "KenCarp47") 
data_cost = data_cost_hessian |> 
  filter(Method_info == "Standard") |> 
  select(Model, Time, solver) |> 
  rename("Time_cost" = "Time")

data_hessian = data_cost_hessian |> 
  filter(Method_info == "ForwardDiff") |> 
  select(Model, Time, solver) |> 
  rename("Time_hessian" = "Time")
data_parameters = data_plot_scaling |> select(Model, n_parameters) |> distinct()
data_plot_hessian = data_hessian |> 
  inner_join(data_cost, by = c("Model", "solver")) |> 
  inner_join(data_parameters, by=c("Model")) |> 
  mutate(n_parameters2 = n_parameters^2) |> 
  filter(!(Model == "Beer_MolBioSystems2014" & solver == "QNDF")) # Crashes

p4 = ggplot(data_plot_hessian, aes(Model, Time_hessian / Time_cost, fill = solver)) + 
  geom_bar(stat="identity", position = position_dodge2(preserve="single", padding = 0.0)) + 
  scale_fill_manual(values = cbPalette[-c(1, 2)], name = "ODE Solver") +
  scale_x_discrete(limits=pos, name="") + 
  geom_text(aes(y = n_parameters2, label=sprintf("%d", n_parameters2))) + 
  coord_flip() + 
  scale_y_log10(expand=c(0, 0, 0.01, 0.1), name = "") + 
  my_minimal + 
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = "none")  

ratio <- data_plot_hessian$Time_hessian / data_plot_hessian$Time_cost
min_val <- min(ratio) * 0.9
max_val <- max(ratio[!is.infinite(ratio)]) * 1.4
data_line <- tibble(x = seq(min_val, max_val, length.out=100), y = x)
p4_alt <- ggplot(filter(data_plot_hessian, !is.infinite(Time_hessian)), aes(n_parameters2, Time_hessian / Time_cost)) + 
  geom_point(aes(color = solver), size=4.0) + 
  geom_line(data=data_line, mapping=aes(x, y), color = cbPalette[1], linewidth=1.0) + 
  labs(x = "Number of parameters to take hessian on squared", y = "Run time Hessian / Run time cost") + 
  scale_x_log10(expand=c(0, 0), limits=c(min_val, max_val)) + 
  scale_y_log10(expand=c(0, 0), limits=c(min_val, max_val)) + 
  scale_color_manual(values = cbPalette[-1]) + 
  theme_classic(base_size = 18) + 
  theme(panel.grid.major = element_line(linewidth=0.05, color="grey10"))

data_sum <- data_plot_hessian |> 
  group_by(solver) |> 
  summarise(is_smaller = sum(Time_hessian / Time_cost < n_parameters2), 
            ratio = sum(Time_hessian / Time_cost < n_parameters2) / n() * 100)

data_table1 <- data_plot_cost |> 
  select(-absTol, -relTol, -Cost_amici) |> 
  rename("Time cost [s]" = "Time")
data_table2 <- data_plot_gradient |> 
  select(-absTol, -relTol, -Gradient_amici) |> 
  rename("Time gradient [s]" = "Time")

table1 <- inner_join(data_table1, data_table2, by = c("Model", "n_parameters", "solver")) |> 
  rename("Number of parameters" = "n_parameters", 
         "Solver" = "solver") |> 
  mutate(Ratio = `Time gradient [s]` / `Time cost [s]`) |> 
  select(Model, Solver, `Time cost [s]`, `Time gradient [s]`, Ratio, everything())
  
table2 <- data_plot_hessian |> 
  select(-n_parameters2) |> 
  mutate(Ratio = Time_hessian / Time_cost) |> 
  rename("Time Hessian [s]" = "Time_hessian", 
         "Time cost [s]" = "Time_cost", 
         "Solver" = "solver", 
         "Number of parameters"="n_parameters") |> 
  select(Model, Solver, `Time Hessian [s]`, `Time cost [s]`, Ratio, everything())

library(gt)

gt_tab1 <- table1 |> 
  gt(rowname_col = "Model") |> 
  tab_header(title = "Cost and gradient benchmark")
for(model in unique(table1$Model)){
  gt_tab1 <- gt_tab1 |> 
    tab_row_group(label = model, 
                  rows = matches(model)) 
}
gt_tab1 <- gt_tab1 |> 
  fmt(columns = starts_with("Time"), fns = function(x) sprintf("%.2e", x)) |> 
  fmt(columns = starts_with("Model"), fns = function(x) sprintf("")) |> 
  fmt(columns = starts_with("Ratio"), fns = function(x) sprintf("%.2e", x)) 

gt_tab2 <- table2 |> 
  gt(rowname_col = "Model")
for(model in unique(table2$Model)){
  gt_tab2 <- gt_tab2 |> 
    tab_row_group(label = model, 
                  rows = matches(model)) 
}
gt_tab2 <- gt_tab2 |> 
  tab_header(title = "Cost and Hessian benchmark") |> 
  fmt(columns = starts_with("Time"), fns = function(x) sprintf("%.2e", x)) |> 
  fmt(columns = starts_with("Model"), fns = function(x) sprintf("")) |> 
  fmt(columns = starts_with("Ratio"), fns = function(x) sprintf("%.2e", x)) 

dir_save <- "Results/Cost_gradient_hessian/"
if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
ggsave(str_c(dir_save, "Cost.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT*1.4)
ggsave(str_c(dir_save, "Gradient.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT*1.4)
ggsave(str_c(dir_save, "Gradient_forward.svg"), p2_, width = BASE_WIDTH, height = BASE_HEIGHT*1.4)

gtsave(gt_tab1, str_c(dir_save, "Table_cost_gradient.tex"))
gtsave(gt_tab2, str_c(dir_save, "Table_cost_hessian.tex"))

ggsave(str_c(dir_save, "Cost_alt.svg"), p1_alt, width = BASE_WIDTH*1.0, height = BASE_HEIGHT*1.0)
ggsave(str_c(dir_save, "Gradient_alt.svg"), p2_alt, width = BASE_WIDTH*1.0, height = BASE_HEIGHT*1.0)
ggsave(str_c(dir_save, "Scaling_alt.svg"), p3_alt, width = BASE_WIDTH*1.0, height = BASE_HEIGHT*1.0)

ggsave(str_c(dir_save, "Scaling.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT*1.4)
ggsave(str_c(dir_save, "Scaling_hessian.svg"), p4, width = BASE_WIDTH, height = BASE_HEIGHT*1.2)
ggsave(str_c(dir_save, "Scaling_hessian_alt.svg"), p4_alt, width = BASE_WIDTH, height = BASE_HEIGHT*1.0)

# ---------------------------------------------------------------------------------------------------------------
# ForwardDiff chunking and keeping parameters fixed 
# ---------------------------------------------------------------------------------------------------------------
dir_result <- "Master-Thesis/Intermediate/Benchmarks/Cost_grad_hess/"
data  <- read_csv(str_c(dir_result, "Fix_parameters.csv"), col_types = cols()) |> 
  filter(!is.infinite(Time))  

color_use = c("#fc8d59", "#d7301f")

data_Lucarelli <- data |> 
  filter(Model == "Fewer_parammodel_Lucarelli_CellSystems2018") |> 
  mutate(N_param_grad = max(N_param_fixed) + 1 - N_param_fixed) 

p1 <- ggplot(data_Lucarelli, aes(N_param_grad, Time, color = chunk_size)) + 
  geom_point(size=3.0) + 
  geom_smooth() + 
  scale_x_continuous(breaks = seq(from = 0, to = max(data_Lucarelli$N_param_grad), by = 5)) + 
  labs(x = "Number of parameters gradient is taken on", 
       y = "Run time [s]", 
       title = "Lucarelli - ForwardDiff chunking reduces run time") + 
  scale_color_manual(values = color_use, name = "Chunk size") + 
  my_minimal + 
  theme(panel.grid.minor = element_blank(), 
        legend.position = "none")

data_bachmann <- data |> 
  filter(Model == "Fewer_parammodel_Bachmann_MSB2011") |> 
  mutate(N_param_grad = max(N_param_fixed) + 1 - N_param_fixed) 

p2 <- ggplot(data_bachmann, aes(N_param_grad, Time, color = chunk_size)) + 
  geom_point(size=3.0) + 
  geom_smooth() + 
  scale_x_continuous(breaks = seq(from = 0, to = max(data_bachmann$N_param_grad), by = 5)) + 
  labs(x = "Number of parameters gradient is taken on", 
       y = "Run time [s]", 
       title = "Bachmann - ForwardDiff chunking reduces run time") + 
  scale_color_manual(values = color_use, name = "Chunk size") + 
  my_minimal + 
  theme(panel.grid.minor = element_blank(), 
        legend.position = "none")

data_isensee <- data |> 
  filter(Model == "Fewer_parammodel_Isensee_JCB2018") |> 
  mutate(N_param_grad = max(N_param_fixed) + 1 - N_param_fixed) 

p3 <- ggplot(data_isensee, aes(N_param_grad, Time, color = chunk_size)) + 
  geom_point(size=3.0) + 
  geom_smooth() + 
  scale_x_continuous(breaks = seq(from = 0, to = max(data_bachmann$N_param_grad), by = 5)) + 
  labs(x = "Number of parameters gradient is taken on", 
       y = "Run time [s]", 
       title = "Isensee - ForwardDiff chunking reduces run time") + 
  scale_color_manual(values = color_use, name = "Chunk size") + 
  my_minimal + 
  theme(panel.grid.minor = element_blank(), 
        legend.position = "none")

ggsave(str_c(dir_save, "Lucarelli_fix.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT*1.0)
ggsave(str_c(dir_save, "Bachmann_fix.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT*1.4)
ggsave(str_c(dir_save, "Isensee_fix.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT*1.4)

# ---------------------------------------------------------------------------------------------------------------
# ForwardDiff chunking with random parameters 
# ---------------------------------------------------------------------------------------------------------------
data <- read_csv(str_c(dir_result, "Test_chunks.csv"), col_types = cols()) |> 
  filter(!is.infinite(Time))  

## Process at nominal values (parameter indices 1)
# Bachmann
data_bachmann = data |> 
  filter(Model == "model_Bachmann_MSB2011") |> 
  filter(I_parameter == 1) |> 
  group_by(Model, chunk_size) |> 
  summarise(median_time = median(Time))
default_data = data_bachmann |> filter(chunk_size == "Default")
data_plot = data_bachmann |> filter(chunk_size != "Default") |> mutate(chunk_size = as.integer(chunk_size))

p1 = ggplot(data_plot, aes(chunk_size, median_time)) + 
  geom_hline(yintercept = default_data$median_time, linewidth=1.5, color=color_use[2]) + 
  geom_point(size=3.0) + 
  geom_smooth() + 
  labs(x = "Chunk size", y = "Run time [s]", title = "Bachmann model") + 
  my_minimal

# Iseense
data_Iseense = data |> 
  filter(Model == "model_Isensee_JCB2018") |> 
  filter(I_parameter == 1) |> 
  group_by(Model, chunk_size) |> 
  summarise(median_time = median(Time))
default_data = data_Iseense |> filter(chunk_size == "Default")
data_plot = data_Iseense |> filter(chunk_size != "Default") |> mutate(chunk_size = as.integer(chunk_size))

p2 = ggplot(data_plot, aes(chunk_size, median_time)) + 
  geom_hline(yintercept = default_data$median_time, linewidth=1.5, color=color_use[2]) + 
  geom_point(size=3.0) + 
  geom_smooth() + 
  labs(x = "Chunk size", y = "Run time [s]", title = "Iseense model") + 
  my_minimal

# Lucarelli
data_Lucarelli = data |> 
  filter(Model == "model_Lucarelli_CellSystems2018") |> 
  filter(I_parameter == 1) |> 
  group_by(Model, chunk_size) |> 
  summarise(median_time = median(Time))
default_data = data_Lucarelli |> filter(chunk_size == "Default")
data_plot = data_Lucarelli |> filter(chunk_size != "Default") |> mutate(chunk_size = as.integer(chunk_size))

p3 = ggplot(data_plot, aes(chunk_size, median_time)) + 
  geom_hline(yintercept = default_data$median_time, linewidth=1.5, color=color_use[2]) + 
  geom_point(size=3.0) + 
  geom_smooth() + 
  labs(x = "Number of chunks used by ForwardDiff.jl", y = "Run time [s]", title = "Lucarelli model") + 
  my_minimal

ggsave(str_c(dir_save, "Bachmann_exhaustive.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT*1.4)
ggsave(str_c(dir_save, "Iseense_exhaustive.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT*1.4)
ggsave(str_c(dir_save, "Lucarelli_exhaustive.svg"), p3, width = BASE_WIDTH, height = BASE_HEIGHT*1.0)

# For each parameter index and one model get the ranking of different chunk-sizes 
get_rank_chunk_size = function(data, model_use)
{
  data_tmp = data |> 
    filter(Model == model_use) |> 
    filter(chunk_size != "Default") |> 
    mutate(chunk_size = as.integer(chunk_size)) |> 
  group_by(chunk_size, Model, I_parameter, solver) |> 
  summarise(Time = median(Time))
  
  
  data_ret = tibble()
  parameter_indices = unique(data_tmp$I_parameter)
  for(i in 1:length(parameter_indices)){
    data_tmp_tmp = data_tmp |> 
      filter(I_parameter == parameter_indices[i])
    data_tmp_tmp = data_tmp_tmp[order(data_tmp_tmp$Time), ] 
    data_tmp_tmp$Rank = 1:nrow(data_tmp_tmp)
    data_ret = bind_rows(data_ret, data_tmp_tmp)
  }
  return(data_ret)  
}


## Lucarelli
data_Lucarelli = get_rank_chunk_size(data, "model_Lucarelli_CellSystems2018")
ggplot(data_Lucarelli, aes(I_parameter, chunk_size, fill=Rank)) + 
  geom_tile(color="white") + 
  scale_fill_viridis_c(direction=-1, option="C") + 
  scale_y_continuous(expand=c(0, 0), breaks=seq(1, max(data_Lucarelli$chunk_size), by = 1)) + 
  scale_x_continuous(expand=c(0, 0), breaks=seq(1, 30, by = 1)) + 
  labs(x = "Parameter index", y = "Chunk size", title = "Lucarelli model - Ranking different chunk-sizes for random parmeters") + 
  my_theme

data_color = data_Lucarelli |> filter(chunk_size %in% c(11, 43, 62)) |> mutate(chunk_size = as.factor(chunk_size))
ggplot(data_color, aes(I_parameter)) + 
  geom_line(data=data_Lucarelli, mapping=aes(I_parameter, y = Time, group = chunk_size), color=cbPalette[1], linewidth=0.8) + 
  geom_line(aes(y = Time, color = chunk_size), linewidth=2.5) + 
  scale_color_manual(values = cbPalette[-1]) + 
  scale_x_continuous(breaks=seq(1, 30, by = 1)) + 
  labs(x = "Parameter index", y = "Run time [s]", title = "Lucarelli model - Run time different chunk-sizes for random parameters") + 
  my_theme


# Bachmann 
data_Bachmann = get_rank_chunk_size(data, "model_Bachmann_MSB2011")
ggplot(data_Bachmann, aes(I_parameter, chunk_size, fill=Rank)) + 
  geom_tile(color="white") + 
  scale_fill_viridis_c(direction=-1, option="C") + 
  scale_y_continuous(expand=c(0, 0), breaks=seq(1, max(data_Bachmann$chunk_size), by = 1)) + 
  scale_x_continuous(expand=c(0, 0), breaks=seq(1, 30, by = 1)) + 
  labs(x = "Parameter index", y = "Chunk size", title = "Bachmann model - Ranking different chunk-sizes for random parmeters") + 
  my_theme

data_color = data_Bachmann |> filter(chunk_size %in% c(7, 13, 14)) |> mutate(chunk_size = as.factor(chunk_size))
ggplot(data_color, aes(I_parameter)) + 
  geom_line(data=data_Bachmann, mapping=aes(I_parameter, y = Time, group = chunk_size), color=cbPalette[1], linewidth=0.8) + 
  geom_line(aes(y = Time, color = chunk_size), linewidth=2.5) + 
  scale_color_manual(values = cbPalette[-1]) + 
  scale_x_continuous(breaks=seq(1, 30, by = 1)) + 
  labs(x = "Parameter index", y = "Run time [s]", title = "Bachmann model - Run time different chunk-sizes for random parameters") + 
  my_theme


## Isensee
data_Iseense = get_rank_chunk_size(data, "model_Isensee_JCB2018")
ggplot(data_Iseense, aes(I_parameter, chunk_size, fill=Rank)) + 
  geom_tile(color="white") + 
  scale_fill_viridis_c(direction=-1, option="C") + 
  scale_y_continuous(expand=c(0, 0), breaks=seq(1, max(data_Iseense$chunk_size), by = 1)) + 
  scale_x_continuous(expand=c(0, 0), breaks=seq(1, 30, by = 1)) + 
  labs(x = "Parameter index", y = "Chunk size", title = "Bachmann model - Ranking different chunk-sizes for random parmeters") + 
  my_theme

data_color = data_Iseense |> filter(chunk_size %in% c(11, 21, 22)) |> mutate(chunk_size = as.factor(chunk_size))
ggplot(data_color, aes(I_parameter)) + 
  geom_line(data=data_Iseense, mapping=aes(I_parameter, y = Time, group = chunk_size), color=cbPalette[1], linewidth=0.8) + 
  geom_line(aes(y = Time, color = chunk_size), linewidth=2.5) + 
  scale_color_manual(values = cbPalette[-1]) + 
  scale_x_continuous(breaks=seq(1, 30, by = 1)) + 
  labs(x = "Parameter index", y = "Run time [s]", title = "Bachmann model - Run time different chunk-sizes for random parameters") + 
  my_theme

# ---------------------------------------------------------------------------------------------------------------
# Adjoint sensitivity analysis 
# ---------------------------------------------------------------------------------------------------------------
data_julia_adjoint <- read_csv(str_c(dir_result, "Test_adjoint_random.csv"), col_types = cols()) |> 
  group_by(Method_info, Model, solver, I_parameter) |> 
  summarise(run_time = median(Time, na.rm = T)) |> 
  mutate(Model = case_when(str_sub(Model, end=6) == "model_" ~ str_sub(Model, start=7), 
                           T ~ Model)) |> 
  mutate(Method_info = str_replace(Method_info, "autojacvec=", ""))

data_julia_adjoint <- read_csv(str_c(dir_result, "Test_adjoint_random.csv"), col_types = cols()) |> 
  group_by(Method_info, Model, solver, I_parameter) |> 
  summarise(run_time = median(Time, na.rm = T)) |> 
  mutate(Model = case_when(str_sub(Model, end=6) == "model_" ~ str_sub(Model, start=7), 
                           T ~ Model)) |> 
  mutate(Method_info = str_replace(Method_info, "autojacvec=", ""))


data_pypesto_ <- data_pypesto_grad_adj |> 
  mutate(Method_info = "AMICI") |> 
  rename("run_time" = "Gradient_amici_adj") |> 
  select(-absTol, -relTol) |> 
  filter(Model %in% c("Bachmann_MSB2011", "Boehm_JProteomeRes2014", "Lucarelli_CellSystems2018"))
n_param_data <- data_pypesto_ |> 
  select(Model, n_parameters)
run_time_norm <- data_pypesto_ |> 
  select(Model, run_time) |> 
  rename("run_time_amici" = "run_time")


# Load adjoint data from AMICI 
data_boehm <- read_csv(str_c("Intermediate/Parameters_test_gradient/Boehm_JProteomeRes2014_results.csv"), col_types = cols()) |> 
  mutate(model = "Boehm_JProteomeRes2014") |> 
  mutate(i_parameter = 1:length(model))
colnames(data_boehm) <- c("tmp1", "run_time", "tmp2", "failed", "Model", "I_parameter")
data_boehm <- data_boehm |> 
  select(-tmp1, -tmp2) |> 
  mutate("Method_info" = "AMICI", 
         "solver" = "AMICI")

data_Lucarelli <- read_csv(str_c("Intermediate/Parameters_test_gradient/Lucarelli_CellSystems2018_results.csv"), col_types = cols()) |> 
  mutate(model = "Lucarelli_CellSystems2018") |> 
  mutate(i_parameter = 1:length(model))
colnames(data_Lucarelli) <- c("tmp1", "run_time", "tmp2", "failed", "Model", "I_parameter")
data_Lucarelli <- data_Lucarelli |> 
  select(-tmp1, -tmp2) |> 
  mutate("Method_info" = "AMICI", 
         "solver" = "AMICI")

data_Bachmann <- read_csv(str_c("Intermediate/Parameters_test_gradient/Bachmann_MSB2011_results.csv"), col_types = cols()) |> 
  mutate(model = "Bachmann_MSB2011") |> 
  mutate(i_parameter = 1:length(model))
colnames(data_Bachmann) <- c("tmp1", "run_time", "tmp2", "failed", "Model", "I_parameter")
data_Bachmann <- data_Bachmann |> 
  select(-tmp1, -tmp2) |> 
  mutate("Method_info" = "AMICI", 
         "solver" = "AMICI")

data_Smith <- read_csv(str_c("Intermediate/Parameters_test_gradient/Smith_BMCSystBiol2013_results.csv"), col_types = cols()) |> 
  mutate(model = "Smith_BMCSystBiol2013") |> 
  mutate(i_parameter = 1:length(model))
colnames(data_Smith) <- c("tmp1", "run_time", "tmp2", "failed", "Model", "I_parameter")
data_Smith <- data_Smith |> 
  select(-tmp1, -tmp2) |> 
  mutate("Method_info" = "AMICI", 
         "solver" = "AMICI")

data_amici <- bind_rows(data_boehm, data_Lucarelli, data_Bachmann, data_Smith) |> 
  filter(I_parameter < 51)

data_adjoint_julia__ <- data_julia_adjoint |> 
  filter(!Method_info %in% c("ForwardDiff", "FiniteDifferences")) |> 
  mutate(failed = ifelse(is.infinite(run_time), T, F)) |> 
  filter(solver != "QNDF") |> 
  filter(Model != "Isensee_JCB2018")

data_plot <- bind_rows(data_amici, data_adjoint_julia__)
data_fail <- data_plot |> 
  group_by(Model, solver, Method_info) |> 
  summarise(n_fail = sum(failed),
            y_coord = 1.1*max(run_time[!is.infinite(run_time)])) 
data_fail_sum <- data_fail |> 
  group_by(Method_info) |> 
  summarise(rate = median(n_fail / 50 * 100), 
            rate_avg = mean(n_fail / 50 * 100),
            median_fail = median(n_fail), 
            total_fail = sum(n_fail))

data_plot <- inner_join(data_plot, data_fail) |> 
  mutate(Method_info = case_when(Method_info == "InterpolatingAdjoint(EnzymeVJP())" ~ "Interpolation - EnzymeVJP", 
                                 Method_info == "InterpolatingAdjoint(ReverseDiffVJP(false))" ~ "Interpolation - ReverseDiffVJP(false)",
                                 Method_info == "InterpolatingAdjoint(ReverseDiffVJP(true))" ~ "Interpolation - ReverseDiffVJP(true)",
                                 Method_info == "QuadratureAdjoint(EnzymeVJP())" ~ "Quadrature - EnzymeVJP", 
                                 Method_info == "QuadratureAdjoint(ReverseDiffVJP(false))" ~ "Quadrature - ReverseDiffVJP(false)",
                                 Method_info == "QuadratureAdjoint(ReverseDiffVJP(true))" ~ "Quadrature - ReverseDiffVJP(true)",
                                 T ~ Method_info)) |> 
    filter(!is.na(Model))
  
data_fail <- data_fail |>   
  mutate(Method_info = case_when(Method_info == "InterpolatingAdjoint(EnzymeVJP())" ~ "Interpolation - EnzymeVJP", 
                                 Method_info == "InterpolatingAdjoint(ReverseDiffVJP(false))" ~ "Interpolation - ReverseDiffVJP(false)",
                                 Method_info == "InterpolatingAdjoint(ReverseDiffVJP(true))" ~ "Interpolation - ReverseDiffVJP(true)",
                                 Method_info == "QuadratureAdjoint(EnzymeVJP())" ~ "Quadrature - EnzymeVJP", 
                                 Method_info == "QuadratureAdjoint(ReverseDiffVJP(false))" ~ "Quadrature - ReverseDiffVJP(false)",
                                 Method_info == "QuadratureAdjoint(ReverseDiffVJP(true))" ~ "Quadrature - ReverseDiffVJP(true)",
                                 T ~ Method_info)) |> 
  filter(!is.na(Model)) 

data_plot$Model <- factor(data_plot$Model, levels = c("Boehm_JProteomeRes2014", "Bachmann_MSB2011", "Lucarelli_CellSystems2018", "Smith_BMCSystBiol2013"))
data_fail$Model <- factor(data_fail$Model, levels = c("Boehm_JProteomeRes2014", "Bachmann_MSB2011", "Lucarelli_CellSystems2018", "Smith_BMCSystBiol2013"))
p_adjoint <- ggplot(filter(data_plot), aes(Method_info, run_time, fill = solver)) + 
  geom_boxplot() + 
  geom_jitter(aes(color=solver), width=0.1, color="black") + 
  facet_wrap(~Model, scale="free_x") + 
  coord_flip() + 
  labs(y = "Run time [s]", x = "") + 
  scale_fill_manual(values = col_highlight[c(5, 1)], name = "Software") + 
  theme_bw(base_size=20) + 
  theme(legend.position = "bottom")

p_fail <- ggplot(data_fail, aes(Method_info, n_fail)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~Model) + 
  coord_flip() + 
  labs(y = "Number of integration failures", x = "", 
       title = "Adjoint sensitivity integration failures", 
       subtitle = "From 50 random parameter vectors") +
  theme_bw(base_size = 22) + 
  scale_y_continuous(expand=c(0, 0), limits=c(0, 50))

ggsave(str_c(dir_save, "Adjoint_random.svg"), p_adjoint, width = BASE_WIDTH*2.2, height = BASE_HEIGHT*2.0)
ggsave(str_c(dir_save, "Adjoint_fail_bar.png"), p_fail, width = BASE_WIDTH*2.2, height = BASE_HEIGHT*2.2, dpi=300)


# Plot failure data 
library(ggalluvial)
data_fail <- data_plot |> mutate(fail = case_when(failed == T ~ "Yes", 
                                                  T ~ "No"))
data_fail1 <- data_fail |>      
    filter(fail == "Yes")
data_fail2 <- data_fail |>      
    filter(fail == "No")
data_fail <- bind_rows(data_fail1, data_fail2)
         
p_adj_fail <- ggplot(data_fail, aes(axis2=Model, 
                      axis1=Method_info)) +
  geom_alluvium(aes(fill=fail), segments=1000, width = 1/12, alpha=.8) +
  geom_stratum(width = 1/12, fill = cbPalette[1], color = "grey") +
  #geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = col_highlight[c(1, 5)]) +
  theme_void() + 
  scale_x_discrete(expand=c(0, 0)) + 
  theme(legend.position = "none")
ggsave(str_c(dir_save, "Adjoint_fail.svg"), p_adj_fail, width = BASE_WIDTH*2.0, height = BASE_HEIGHT*2.0)


# Look closer where Julia does not fail (that frequently)
data_boehm_closer <- data_plot |> 
  filter(Model == "Boehm_JProteomeRes2014") |> 
  filter(Method_info %in% c("Interpolation - EnzymeVJP", "AMICI"))
data_amici__ <- data_boehm_closer |> 
  filter(Method_info == "AMICI") |> 
  rename("run_time_AMICI" = "run_time") |> 
  select(run_time_AMICI, I_parameter)
data_plot_boehm_closer <- data_boehm_closer |> 
  filter(Method_info != "AMICI") |> 
  inner_join(data_amici__)

max_val <- max(data_plot_boehm_closer$run_time_AMICI)*1.05
range_data <- tibble(x = c(0.005, max_val), y = c(0.005, max_val))
p1 <- ggplot(filter(data_plot_boehm_closer, !is.infinite(run_time)), aes(run_time, run_time_AMICI)) + 
  geom_point(size=5.0) +
  geom_abline(slope=1.0, linewidth=2.0, color = cbPalette[1]) + 
  geom_abline(slope=1.5, linetype=2, linewidth=2.0, color = cbPalette[1]) + 
  geom_abline(slope=2.0, linetype=2, linewidth=2.0, color = cbPalette[1]) + 
  xlim(0.005, max_val) + ylim(0.005, max_val) +
  labs(y = "Run time AMICI [s]", x = "Run time PEtab.jl [s]") + 
  theme_classic(base_size = 22) + 
  theme(panel.grid.major.y = element_blank(), legend.position = "none") +
  theme(axis.text = element_text(size=24))


data_bachmann_closer <- data_plot |> 
  filter(Model == "Bachmann_MSB2011") |> 
  filter(Method_info %in% c("Interpolation - EnzymeVJP", "AMICI"))
data_amici__ <- data_bachmann_closer |> 
  filter(Method_info == "AMICI") |> 
  rename("run_time_AMICI" = "run_time") |> 
  select(run_time_AMICI, I_parameter)
data_plot_bachmann_closer <- data_bachmann_closer |> 
  filter(Method_info != "AMICI") |> 
  inner_join(data_amici__)

max_val <- max(data_plot_bachmann_closer$run_time_AMICI[!is.infinite(abs(data_plot_bachmann_closer$run_time_AMICI))])*1.05
range_data <- tibble(x = c(0.1, max_val), y = c(0.1, max_val))
p2 <- ggplot(filter(data_plot_bachmann_closer, !is.infinite(run_time)), aes(run_time, run_time_AMICI)) + 
  geom_point(size=5.0) +
  geom_abline(slope=1.0, linewidth=2.0, color = cbPalette[1]) + 
  geom_abline(slope=1.5, linetype=2, linewidth=2.0, color = cbPalette[1]) + 
  geom_abline(slope=2.0, linetype=2, linewidth=2.0, color = cbPalette[1]) + 
  xlim(0.1, max_val) + ylim(0.1, max_val) +
  labs(y = "Run time AMICI [s]", x = "Run time PEtab.jl [s]") + 
  theme_classic(base_size = 22) + 
  theme(panel.grid.major.y = element_blank(), legend.position = "none") +
  theme(axis.text = element_text(size=24))

val_tmp1 <-  (data_plot_boehm_closer |> 
  filter(!is.infinite(run_time)) |> 
  mutate(ratio = 1 / mean(run_time / run_time_AMICI)))$ratio[1]
val_tmp2 <-  (data_plot_bachmann_closer |> 
  filter(!is.infinite(run_time)) |> 
  mutate(ratio = 1 / mean(run_time / run_time_AMICI)))$ratio[1]

ggsave(str_c(dir_save, "Adjoint_boehm_closer.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(str_c(dir_save, "Adjoint_bachmann_closer.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)

# Check accuracy of gradients Boehm
model_name <- "Boehm_JProteomeRes2014"
data_julia_adjoint <- read_csv(str_c(dir_result, "Test_adjoint_random_gradientmodel_Boehm_JProteomeRes2014.csv"), col_types = cols()) |> 
  group_by(Method_info, Model, solver, I_parameter) |> 
  mutate(Model = case_when(str_sub(Model, end=6) == "model_" ~ str_sub(Model, start=7), 
                           T ~ Model)) |> 
  mutate(Method_info = str_replace(Method_info, "autojacvec=", ""))
data_amici <- read_csv("./Intermediate/Parameters_test_gradient/Boehm_JProteomeRes2014_results_grad.csv", 
                       col_types = cols()) |> 
  mutate(I_parameter = 1:100) |> 
  filter(I_parameter < 51) |> 
  mutate("Method_info" = "AMICI", 
         "solver" = "AMICI", 
         model = model_name)
parameters <- read_csv("./Intermediate/Parameters_test_gradient/Boehm_JProteomeRes2014.csv", 
                       col_types = cols()) 

data_amici <- data_amici[, 2:ncol(data_amici)]
data_julia_adjoint <- bind_rows(data_julia_adjoint, data_amici)

n_param <- 9
reference_gradient = data_julia_adjoint |> 
  filter(Method_info == "ForwardDiff")
methods_grad <- c("InterpolatingAdjoint(ReverseDiffVJP(true))", "InterpolatingAdjoint(ReverseDiffVJP(false))", "InterpolatingAdjoint(EnzymeVJP())", "QuadratureAdjoint(ReverseDiffVJP(true))", "QuadratureAdjoint(ReverseDiffVJP(false))", "QuadratureAdjoint(EnzymeVJP())", "AMICI")
data_grad_accuracy <- tibble()
for(i in 1:50){
  if(i == 1) next
  
  data_ref <- (reference_gradient |> 
                 filter(I_parameter == i) |> 
                 unlist(use.names = F))[1:n_param] |> 
    as.numeric()  
  for(j in 1:length(methods_grad)){
    data_test <- data_julia_adjoint |> 
      filter(I_parameter == i & Method_info == methods_grad[j])
    data_test <- as.numeric((data_test |> unlist(use.names = F))[1:n_param])
    if(sum(data_test) == 0 || is.infinite(sum(data_test))){
      norm_val <- Inf
    }else if(methods_grad[j] == "AMICI"){
      parameter_values <- as.numeric(parameters[i, ] |> unlist(use.names = F))
      data_test__ <- 10^parameter_values * log(10) * data_test * (-1)
      norm_val <- sqrt(sum((data_ref - data_test__)^2))
    }else{
      norm_val <- sqrt(sum((data_ref - data_test)^2))
    }
    data_tmp <- tibble(I_parameter = i, 
                       method_info = methods_grad[j], 
                       model = model_name, 
                       L2_norm = norm_val / n_param)    
    data_grad_accuracy <- bind_rows(data_grad_accuracy, data_tmp)
  }
}
data_grad_accuracy_boehm <- data_grad_accuracy |> 
  mutate(is_amici = case_when(method_info == "AMICI" ~ "Yes", 
                              T ~ "No"))

# Check accuracy of gradients Boehm
model_name <- "Bachmann_MSB2011"
data_julia_adjoint <- read_csv(str_c(dir_result, "Test_adjoint_random_gradientmodel_Bachmann_MSB2011.csv"), col_types = cols()) |> 
  group_by(Method_info, Model, solver, I_parameter) |> 
  mutate(Model = case_when(str_sub(Model, end=6) == "model_" ~ str_sub(Model, start=7), 
                           T ~ Model)) |> 
  mutate(Method_info = str_replace(Method_info, "autojacvec=", ""))
data_amici <- read_csv("./Intermediate/Parameters_test_gradient/Bachmann_MSB2011_results_grad.csv", 
                       col_types = cols()) |> 
  mutate(I_parameter = 1:100) |> 
  filter(I_parameter < 51) |> 
  mutate("Method_info" = "AMICI", 
         "solver" = "AMICI", 
         model = model_name)
parameters <- read_csv("./Intermediate/Parameters_test_gradient/Bachmann_MSB2011.csv", 
                       col_types = cols()) 
data_amici <- data_amici[, 2:ncol(data_amici)]
data_julia_adjoint <- bind_rows(data_julia_adjoint, data_amici)

n_param <- 113
reference_gradient = data_julia_adjoint |> 
  filter(Method_info == "ForwardDiff")
methods_grad <- c("InterpolatingAdjoint(ReverseDiffVJP(true))", "InterpolatingAdjoint(ReverseDiffVJP(false))", "InterpolatingAdjoint(EnzymeVJP())", "QuadratureAdjoint(ReverseDiffVJP(true))", "QuadratureAdjoint(ReverseDiffVJP(false))", "QuadratureAdjoint(EnzymeVJP())", "AMICI")
data_grad_accuracy <- tibble()

for(i in 1:50){
  if(i == 1) next
  
  data_ref <- (reference_gradient |> 
                 filter(I_parameter == i) |> 
                 unlist(use.names = F))[1:n_param] |> 
    as.numeric()  
  for(j in 1:length(methods_grad)){
    data_test <- data_julia_adjoint |> 
      filter(I_parameter == i & Method_info == methods_grad[j])
    data_test <- as.numeric((data_test |> unlist(use.names = F))[1:n_param])
    if(sum(data_test) == 0 || is.infinite(sum(data_test))){
      norm_val <- Inf
    }else if(methods_grad[j] == "AMICI"){
      parameter_values <- as.numeric(parameters[i, ] |> unlist(use.names = F))
      data_test__ <- 10^parameter_values * log(10) * data_test * (-1)
      norm_val <- sqrt(sum((data_ref - data_test__)^2))
    }else{
      norm_val <- sqrt(sum((data_ref - data_test)^2))
    }
    data_tmp <- tibble(I_parameter = i, 
                       method_info = methods_grad[j], 
                       model = model_name, 
                       L2_norm = norm_val / n_param)    
    data_grad_accuracy <- bind_rows(data_grad_accuracy, data_tmp)
  }
}
data_grad_accuracy_bachmann <- data_grad_accuracy |> 
  mutate(is_amici = case_when(method_info == "AMICI" ~ "Yes", 
                              T ~ "No"))

ggplot(filter(data_grad_accuracy_boehm, !is.infinite(L2_norm)), aes(method_info, L2_norm, fill = is_amici)) + 
  geom_violin(draw_quantiles=0.5, linewidth=1.0) +
  geom_jitter(alpha=1.0, width = 0.2) + 
  scale_y_log10() + 
  scale_fill_manual(values = col_highlight[c(1, 5)]) + 
  my_theme

ggplot(filter(data_grad_accuracy_bachmann, !is.infinite(L2_norm)), aes(method_info, L2_norm, fill = is_amici)) + 
  geom_violin(draw_quantiles=0.5, linewidth=1.0) +
  geom_jitter(alpha=1.0, width = 0.2) + 
  scale_y_log10() + 
  scale_fill_manual(values = col_highlight[c(1, 5)]) + 
  my_theme

data_boehm_enzyme <- data_grad_accuracy_boehm |> 
  filter(method_info == "InterpolatingAdjoint(EnzymeVJP())" | method_info == "AMICI")
data_bachmann_enzyme <- data_grad_accuracy_bachmann |> 
  filter(method_info == "InterpolatingAdjoint(EnzymeVJP())" | method_info == "AMICI")

p1 <- ggplot(filter(data_boehm_enzyme, !is.infinite(L2_norm)), aes(method_info, L2_norm)) + 
  geom_boxplot(fill = cbPalette[1]) +
  geom_jitter(alpha=1.0, width = 0.2) +
  labs(x = "", y = "") +
  scale_y_log10() + 
  my_theme

p2 <- ggplot(filter(data_bachmann_enzyme, !is.infinite(L2_norm)), aes(method_info, L2_norm)) + 
  geom_boxplot(fill = cbPalette[1]) +
  geom_jitter(alpha=1.0, width = 0.2) + 
  labs(x = "", y = "") +
  scale_y_log10() + 
  my_theme

ggsave(str_c(dir_save, "Accuracy_boehm_enzyme.svg"), p1, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave(str_c(dir_save, "Accuracy_bachmann_enzyme.svg"), p2, width = BASE_WIDTH, height = BASE_HEIGHT)
