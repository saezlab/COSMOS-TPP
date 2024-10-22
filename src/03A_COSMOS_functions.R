download_omnipath <- function(url){
  
  return(read.table(url, sep = '\t', header = TRUE))
}

format_dorothea <- function(url = 'http://omnipathdb.org/interactions?datasets=tfregulons&tfregulons_levels=A,B,C&genesymbols=1&fields=sources,tfregulons_level')
{

  dorothea <- download_omnipath(url)
  dorothea <- dorothea[,c(4,3,6,7)]
  dorothea$sign <- dorothea$is_stimulation - dorothea$is_inhibition
  dorothea <- dorothea[dorothea$sign != 0,]
  dorothea <- dorothea[,c(1,2,5)]
  
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  G_list <- getBM(filters = "hgnc_symbol",
                  attributes = c('hgnc_symbol','entrezgene_id', "description"),
                  values = unique(dorothea$source_genesymbol), mart = ensembl)
  
  gene_mapping <- as.character(G_list[,2])
  names(gene_mapping) <- G_list[,1]
  
  for(i in 1:length(dorothea[,1]))
  {
    dorothea[i,2] <- gene_mapping[dorothea[i,2]]
  }
  
  dorothea$source_genesymbol <- paste0("X",dorothea$source_genesymbol)
  
  return(dorothea)
}

format_dorothea_pck <- function()
{
  
  dorothea <- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A","B","C","D"),c(3,1,4)])
  names(dorothea) <- c("target_genesymbol","source_genesymbol","sign")
  
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  
  G_list <- getBM(filters = "hgnc_symbol",
                  attributes = c('hgnc_symbol','entrezgene_id', "description"),
                  values = unique(dorothea$source_genesymbol), mart = ensembl)
  
  gene_mapping <- as.character(G_list[,2])
  names(gene_mapping) <- G_list[,1]
  
  for(i in 1:length(dorothea[,1]))
  {
    dorothea[i,2] <- gene_mapping[dorothea[i,2]]
  }
  
  dorothea$source_genesymbol <- paste0("X",dorothea$source_genesymbol)
  
  return(dorothea)
}

####################################################################################################
####################################################################################################
####################################################################################################

downstream_neighbours <- function(meta_network, n_steps, input_names)
{
  meta_g <- graph_from_data_frame(meta_network[,c(1,3,2)]) 
  
  dn_nbours <- ego(graph = meta_g, order = n_steps,nodes = input_names,mode = "out")
  
  sub_nodes <- unique(names(unlist(dn_nbours)))
  
  meta_network <- meta_network[meta_network$source %in% sub_nodes & meta_network$target %in% sub_nodes,]
  
  return(meta_network)
}

filter_inputs <- function(inputs, meta_network)
{
  meta_g <- graph_from_data_frame(meta_network[,c(1,3,2)]) 
  
  inputs <- inputs[,names(inputs) %in% V(meta_g)$name]
  
  return(inputs)
}

####################################################################################################
####################################################################################################
####################################################################################################

filter_TF_sign <- function(meta_network, ttop_RNA, inputs, TF_targets, my_threshold = 1)
{
  TF_targets <- TF_targets[TF_targets$source_genesymbol %in% meta_network$source,]
  
  meta_network$target_sign <- ifelse(meta_network$target %in% ttop_RNA$ID, ttop_RNA$FC, NA) 
  
  inputs <- inputs[,names(inputs) %in% TF_targets[,2]]
  
  meta_network$TF_sign <- sapply(meta_network$source, function(x, inputs)
  {
    if(x %in% names(inputs))
    {
      return(as.numeric(inputs[,x]))
    } else
    {
      return(NA)
    }
  },inputs = inputs)
  
  meta_network$target_sign <- ifelse(abs(meta_network$target_sign) > my_threshold, meta_network$target_sign, 0)
  
  save_df <- meta_network[(which(is.na(!(meta_network$source %in% TF_targets[,2] & meta_network$target_sign == 0)))),]
  
  meta_network <- meta_network[!(meta_network$source %in% TF_targets[,2] & meta_network$target_sign == 0),]
  
  meta_network <- meta_network[is.na(meta_network$TF_sign) | sign(meta_network$target_sign*meta_network$interaction) == sign(meta_network$TF_sign),]
  meta_network <- meta_network[,c(1,2,3)]
  
  meta_network <- meta_network[complete.cases(meta_network),]
  meta_network <- as.data.frame(rbind(meta_network,save_df[,c(1,2,3)]))
  
  return(meta_network)
}

get_TF_sign_from_CARNI <- function(carni_res, TF_targets)
{
  TF_signs <- as.data.frame(carni_res$nodesAttributes)
  TF_signs <- TF_signs[TF_signs$AvgAct != 0,]
  
  TF_signs <- TF_signs[TF_signs$Node %in% TF_targets[,2],]
  row.names(TF_signs) <- TF_signs$Node
  TF_signs <- TF_signs[,"AvgAct",drop = F]
  TF_signs <- as.data.frame(t(TF_signs))
  
  return(TF_signs)
}

makeSubNet <- function(source, mode, sif, att)
{
  ig_net <- graph_from_data_frame(sif) 
  subnet <- shortest_paths(ig_net, from = source, to = att[att$type == "metabolite",1], mode = mode)
  
  subnet <- sif[sif$Node1 %in% names(unlist(subnet$vpath)) & sif$Node2 %in% names(unlist(subnet$vpath)),]
  
  att <- att[att$Nodes %in% subnet$Node1 | att$Nodes %in% subnet$Node2,]
  
}

display_node_neighboorhood <- function(central_node,sif, att, n = 100)
{
  full_sif <- sif
  full_att <- att
  
  ig_net <- graph_from_data_frame(full_sif) 
  
  ig_net <- make_ego_graph(ig_net, nodes = central_node, order = n, mode = "all")[[1]]
  
  to_keep <- V(ig_net)$name
  
  full_sif <- full_sif[full_sif$Node1 %in% to_keep & full_sif$Node2 %in% to_keep,]
  full_att <- full_att[full_att$Nodes %in% to_keep,]
  
  center_node <- shortest_paths(ig_net, from = central_node, to = full_att[full_att$measured == 1,1])
  center_node_out <- full_sif[full_sif$Node1 %in% names(unlist(center_node$vpath)) & full_sif$Node2 %in% names(unlist(center_node$vpath)),]
  
  # write_csv(center_node_net,"center_node_sif_newDoro.csv")
  
  center_node <- shortest_paths(ig_net, from = central_node, to = full_att[full_att$measured == 1,1], mode = "in")
  center_node_in <- full_sif[full_sif$Node1 %in% names(unlist(center_node$vpath)) & full_sif$Node2 %in% names(unlist(center_node$vpath)),]
  
  center_node_net <- as.data.frame(rbind(center_node_in,center_node_out))
  
  nodes <- full_att[full_att$Nodes %in% center_node_net$Node1 | full_att$Nodes %in% center_node_net$Node2,]
  edges <- center_node_net
  
  names(edges) <- c("from","to","sign","weigth")
  edges$color <- ifelse(edges$sign == 1, "green","red")
  edges$arrows <- "to"
  edges <- unique(edges)
  
  names(nodes)[1] <- "id"
  nodes$label <- nodes$id
  nodes$color <- ifelse(nodes$Activity > 0, "green","red")
  nodes <- nodes[!duplicated(nodes$id),]
  nodes$shape <- "dot" 
  nodes[nodes$type == "metab_enzyme","shape"] <- "square"
  nodes[nodes$type == "protein","shape"] <- "square"
  nodes[nodes$type == "Kinase","shape"] <- "triangle"
  nodes[nodes$type == "TF","shape"] <- "diamond"
  nodes <- nodes[order(nodes$id),]
  nodes$shadow <- ifelse(nodes$measured == 1, T, F)
  
  return(visNetwork(nodes, edges, width = "100%") %>% 
           visOptions(highlightNearest = TRUE, 
                      nodesIdSelection = list(enabled = TRUE,
                                              style = 'width: 200px; height: 26px;
                                              background: #f8f8f8;
                                              color: darkblue;
                                              border:none;
                                              outline:none;')))
  
}

df_to_viper_regulon <- function(df)
{
  names(df) <- c("feature","pathway","sign")
  df <- df[complete.cases(df),]
  
  pathway_regulon <- list(0)
  i <- 1
  for(pathway in unique(df$pathway))
  {
    pathway_feature_list <- list(0)
    features <- df[df$pathway == pathway, 3]
    names(features) <- df[df$pathway == pathway, 1]
    pathway_feature_list[[1]] <- features
    pathway_feature_list[[2]] <- rep(1,length(features))
    names(pathway_feature_list) <- c("tfmode","likelihood")
    
    pathway_regulon[[i]] <- pathway_feature_list
    i <- i+1
  }
  names(pathway_regulon) <- unique(df$pathway)
  return(pathway_regulon)
}

is_expressed <- function(x)
{
  if(!grepl("Metab",x))
  {
    if(gsub("X","",x) %in% expressed_gene_list)
    {
      return(x)
    } else
    {
      if(grepl("XGene[0-9]+__[0-9_]+$",x))
      {
        genes <- gsub("XGene[0-9]+__","",x)
        genes <- strsplit(genes,"_")[[1]]
        if(sum(genes %in% expressed_gene_list) != length(genes))
        {
          return(NA)
        } else
        {
          return(x)
        }
      } else
      {
        if(grepl("XGene[0-9]+__[0-9_]+reverse",x))
        {
          genes <- gsub("XGene[0-9]+__","",x)
          genes <- gsub("_reverse","",genes)
          genes <- strsplit(genes,"_")[[1]]
          if(sum(genes %in% expressed_gene_list) != length(genes))
          {
            return(NA)
          } else
          {
            return(x)
          }
        } else
        {
          return(NA)
        }
      }
    } 
  } else
  {
    return(x)
  }
}
