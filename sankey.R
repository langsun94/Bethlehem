#source tracker sankey diagram
library(networkD3)

json <- "/Volumes/efs/CANR/Ansci/D\'Amico\ Lab/Langsun_las17015/2.\ Bethlehem\ project/Bethlehem2020/Bacteria/source\ tracker/sankey/bacteria.json"
st <- jsonlite::fromJSON(json)

group <- as.character(c(1,0,0,0,
                        2,0,3,4,
                        4,4,5,6,
                        7,8,8,9,
                        9,10))

sankeyNodes <- data.frame(st$nodes, group)

p = sankeyNetwork(Links = st$links, Nodes = sankeyNodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              units = "TWh", fontSize = 12, nodeWidth = 30,
              sinksRight=F,iterations = 0,NodeGroup="group")

library(htmlwidgets)
saveWidget(p, file=paste0( getwd(), "sankey.html"))


