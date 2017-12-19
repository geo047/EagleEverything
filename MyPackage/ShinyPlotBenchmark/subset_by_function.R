df <- read.csv(file="\\\\braggflush1\\flush1\\bow355\\AMplus_new_code\\Mid_docker_tests\\benchmark\\new_timing.out")
# df <- read.csv(file="\\\\braggflush1\\flush1\\bow355\\AMplus_new_code\\Large\\new_timing.out")
df[ order(df$function., df$ncpu, df$ngpu), ]

df_sub <- subset(df, ngpu==0 & itnum>1 , c(function.,ncpu,time_ms))

df_sub[order(df_sub$ncpu),]

df_sub <- df_sub[order(df_sub$function.,df_sub$ncpu),]

df_sub_calc_extBIC <- subset(df_sub, function.=="calc_extBIC" , c(ncpu,time_ms))


library(rbokeh)

figure() %>%
    ly_points(x=ncpu, y=time_ms, data = df_sub, color = function., hover = list(time_ms,function.,ncpu)) %>%
    y_axis(label = "Time/s", log=F) %>%
    y_range(c(-100, max(df_sub$time_ms)+0.1*df_sub$time_ms)  ) %>%
    x_axis(label = "Eagle Function") 
    
    
#   ly_bar( x= as.factor(df_sub_calc_extBIC$ncpu), y= df_sub_calc_extBIC$time_ms,
#            position = "fill", width = 2) 
