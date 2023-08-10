library(ggplot2)
library(reshape2)
df <- read.csv('/home/mohamed/Data/Work/Prairie_algorithm/2-HDS_v2/HDS_standalone_code/HDS_output.csv') 

df <- df[,-ncol(df)]

df_plt <- melt(data = df, c('time', 'basinID'))
df_plt <- df_plt[df_plt$variable=='volFrac' | df_plt$variable=='conArea' ,]
ggplot(data = df_plt, aes(x=time, y=value))+
  geom_line()+
  facet_grid(basinID ~ variable, scales = 'free_y')


ggplot(df[df$basinID==3,],aes(x = volFrac,y = conArea)) +
  geom_segment(aes(xend=c(tail(volFrac, n=-1), NA), yend=c(tail(conArea, n=-1), NA)),size=0.5,
               arrow=arrow(length=unit(0.05,"in"),type = 'open'))+
  #facet_zoom(zoom.data = hydroyr ==2013,xlim = c(10,20))+ #x = which(data2$ZPND_mm>5 & data2$ZPND_mm<10))+#,ylim = c(0.04,0.06))+
  facet_wrap(basinID~., nrow = 1)+theme(aspect.ratio = 1)
