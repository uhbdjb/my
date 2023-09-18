mydata.markers <- FindAllMarkers(mydata,only.pos = TRUE, win.pct =0.25, logfc.threshold = 1)
write.csv(mydata.markers,file='marker.csv')

a<-mydata.markers %>% group_by(cluster) %>% top_n(n =2,wt = avg_log2FC)

top5 <- mydata.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(mydata, features = top5$gene,label=F,group.bar.height = 0.03, group.colors = c( "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7"))+ scale_fill_gradientn(colors = c( "navy", "white", "firebrick3"))
