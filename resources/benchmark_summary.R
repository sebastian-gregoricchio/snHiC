benchmark_summary = 
  function(bechnmark_dir) {
    require("dplyr")
    require("lubridate")
    
    benchmark_dir = gsub("/$", "", benchmark_dir)
    
    tables = list.files(path = benchmark_dir, recursive = T, full.names = T)
    tb_list = lapply(tables, function(x){data.frame(data.table::fread(x))})
    
    summary_tb =
      data.frame(dplyr::mutate(do.call(rbind, tb_list),
                               process = basename(tables)) %>%
                   tidyr::separate(col = process,
                                   into = c("Rule", "Step"),
                                   sep = "-",
                                   fill = "right") %>%
                   purrr::map_df(~ gsub("[$-]", 0, .x)) %>%
                   dplyr::mutate(Rule = gsub("[.]tsv","",Rule),
                                 Step = gsub("[.]tsv","",Step),
                                 s = as.numeric(s),
                                 max_rss = as.numeric(max_rss),
                                 max_vms = as.numeric(max_vms),
                                 max_uss = as.numeric(max_uss),
                                 max_pss = as.numeric(max_pss),
                                 io_in = as.numeric(io_in),
                                 io_out = as.numeric(io_out),
                                 mean_load = as.numeric(mean_load),
                                 cpu_time = as.numeric(cpu_time)) %>%
                   dplyr::group_by(Rule) %>%
                   dplyr::summarise(N.steps = n(),
                                    Tot_Running_Time_min = round(sum(s, na.rm = T)/60, 1),
                                    Tot_Running_Time_dd.hh.mm.ss = gsub("[.].*$","s",lubridate::seconds_to_period(sum(s, na.rm = T))),
                                    Max_physical_mem_GB = round(max(max_rss, na.rm = T)/1024, 1),
                                    Max_virtual_mem_GB = round(max(max_vms, na.rm = T)/1024, 1),
                                    Average_mean.load = round(mean(mean_load, na.rm = T), 1)) %>%
                   dplyr::arrange(Rule))
    
    total = data.frame(Rule = "SUMMARY",
                       N.steps = sum(summary_tb$N.steps, na.rm = T),
                       Tot_Running_Time_min = round(sum(summary_tb$Tot_Running_Time_min, na.rm = T), 1),
                       Tot_Running_Time_dd.hh.mm.ss = gsub("[.].*$","s",lubridate::seconds_to_period(sum(summary_tb$Tot_Running_Time_min, na.rm = T)*60)),
                       Max_physical_mem_GB = round(max(summary_tb$Max_physical_mem_GB, na.rm = T), 1),
                       Max_virtual_mem_GB = round(max(summary_tb$Max_virtual_mem_GB, na.rm = T), 1),
                       Average_mean.load = round(mean(summary_tb$Average_mean.load, na.rm = T), 1))
    
    return(rbind(summary_tb, total))
  }
