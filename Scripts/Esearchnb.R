library(Esearch3D)


#1) A two column dataframe representing the gene - fragment interaction network
#gf_nb <- gf_nb %>% semi_join(ann_net_nb, by = c("V2" = "ID"))
gf_net = as.matrix(gf_nb)
#2) A two column dataframe representing the fragment-fragment interaction network
ff_net = as.matrix(nBff)

#3) A matrix containing the starting values of the nodes (gene node scores are given by the RNA-seq values while fragments have a starting value of 0)
input_m = as.matrix(inputnB)
#4) A list of genes to annotate to the chromatin fragments
ann_net_b= ann_net_nb

#Two step propagation -----
#Propagation over the gene-fragment network
gf_prop = rwr_OVprop(g=gf_net,
                     input_m = input_m, 
                     no_cores = 5, 
                     r = 0.8)
#Propagation over the gene-fragment-fragment network
ff_prop=rwr_OVprop(g=ff_net,
                   input_m = gf_prop, 
                   no_cores = 5, 
                   r = 0.8)

#Create an annotated igraph object
net=create_net2plot(gf_net,
                    input_m,
                    gf_prop,
                    ann_net_b,
                    frag_pattern = "frag",
                    ff_net,
                    ff_prop)

#Start the GUI
start_GUI(net, ann_net_b)
      
