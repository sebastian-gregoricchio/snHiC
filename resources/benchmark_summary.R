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
                   dplyr::mutate(Rule = gsub("[.]tsv","",Rule),
                                 Step = gsub("[.]tsv","",Step)) %>%
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
                       Tot_Running_Time_min = sum(summary_tb$Tot_Running_Time_min, na.rm = T),
                       Tot_Running_Time_dd.hh.mm.ss = gsub("[.].*$","s",lubridate::seconds_to_period(sum(summary_tb$Tot_Running_Time_min, na.rm = T)*60)),
                       Max_physical_mem_GB = max(summary_tb$Max_physical_mem_GB, na.rm = T),
                       Max_virtual_mem_GB = max(summary_tb$Max_virtual_mem_GB, na.rm = T),
                       Average_mean.load = mean(summary_tb$Average_mean.load, na.rm = T))
    
    return(rbind(summary_tb, total))
  }
