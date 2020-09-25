library(KEGGgraph) # trans_xml2graph
library(graphite) # get pathway
library(igraph)
library(xlsx)
library(BioNet) # get edgelist



# download kgml files from https://www.kegg.jp/kegg/xml/
# translate the kgml files to graph object 

# spec is the species name in kegg ("hsa","mmu",etc.)
# file_root is the root you want to save kegg pathway kgml files

 
# download kegg kgml files
kegg_download <- function(spec,file_root){
  # create the file dir for downloading
  if (!dir.exists(paste0(file_root,"/kgml"))){
    dir.create(paste0(file_root,"/kgml"))
  }
  download_dir <- paste0(file_root,"/kgml/",spec)
  if (!dir.exists(download_dir)){
    dir.create(download_dir)
  }

  url_pathway_names <- paste0("http://rest.kegg.jp/list/pathway/",spec)
  # downloaded pathway info files
  download.file(url = url_pathway_names, destfile = paste0(file_root,"/",spec,"_pathways.txt"))
  file1 = read.csv(file = paste0(file_root,"/",spec,"_pathways.txt"),header = F,sep = "\t")
  if (length(unique(substr(as.vector(file1$V1),1,5)))!=1){
    stop("You downloaded the wrong pathway info file! Please check your network connection!")
  }
  name <- substr(as.vector(file1$V1),6,13)
  # repeat several times for complete download
  download_repeat <- function(name,download_dir){
    loaded_files <- list.files(download_dir)
    name_undownload <- name[!sapply(X = name,FUN = function(x)(paste0(x,".xml") %in% loaded_files))]
    error_files <- sapply(name_undownload,function(x)("try-error" %in% class(download.file(url = paste0("http://rest.kegg.jp/get/",x,"/kgml"),
                                                                                           destfile = paste0(download_dir,"/",x,".xml"),
                                                                                           quiet=TRUE))))
    if (length(name_undownload)!=0){
      name_undownload <- name_undownload[error_files]
    }
    name_undownload
  }
  
  print(paste0("Downloading ",spec," KEGG pathays..."))
  print("This may take you ten minutes, depending on your internet speed.")
  error_file <- download_repeat(name,download_dir)
  # try 5 times for completely downloading
  repeat_num = 5
  while (length(error_file)!=0 | repeat_num == 0){
    error_file <- download_repeat(name,download_dir)
    repeat_num <- repeat_num-1
  }
  if (length(error_file)!=0){
    stop("Incomplete download! Please try again later or manually download.",
         "\n These pathways failed to download:",error_file)
  }else{
    print(paste0("You have completely downloaded all ",spec," kgml files in kegg pathways!"))
  }
}


# translate kgml files to graphNEl
trans_xml2graph <- function(xml_list,file_dir,save_dir="."){
  graph_list <- sapply(X = xml_list, FUN = function(file_name) try(parseKGML2Graph(paste0(file_dir,"/",file_name))))
  # remove pathways which failed to be translated
  right_trans <- sapply(graph_list,function(x)is(x,"graphNEL"))
  wrong_trans_file <- names(graph_list)[!right_trans]
  if (length(wrong_trans_file)!=0){
    print("These files failed to translate!")
    print(wrong_trans_file)
  }
  graph_list <- graph_list[right_trans]
  # remove pathways which have no elements in nodes or edges
  zero_id <- sapply(X = graph_list,FUN = function(x)(length(nodes(x))==0 | length(edgeData(x))==0))
  graph_list <- graph_list[!zero_id]
  names(graph_list) <- t(data.frame(strsplit(x = names(graph_list),split = "\\.")))[,1]
  # remove the specices label in ENTREZID  gene id
  fun_trans <- function(id,graph_list){
    graph0 <- graph_list[[id]]
    nodes(graph0) <- as.vector(sapply(nodes(graph0), function(x)gsub(pattern = paste0(spec,":"),replacement = "",x = x)))
    return(graph0)
  }
  graph_list <- sapply(names(graph_list), FUN = fun_trans,graph_list)
  save(list = c("graph_list"),file = paste0(save_dir,"/graph_list.RData"))
  graph_list
}


# get kegg pathways' edges and nodes information
get_kegg_info <- function(graph_list, save_dir=NULL){
  kle <- lapply(graph_list, getEdgeList)
  kleu <- lapply(kle,function(x)x[c("from","to")])
  kegg_edges <- lapply(kleu, function(x)x[!duplicated(x),])
  # kegg_edges <- do.call('rbind',kegg_edges)
  kegg_nodes <- sapply(graph_list,function(x)nodes(x))
  if (!is.null(save_dir)){
    if (substr(save_dir,nchar(save_dir),nchar(save_dir)) == "/"){
      save_dir0 <- substr(save_dir,1,(nchar(save_dir)-1))
    }else{
      save_dir0 <- save_dir
    }
    save(list = c("kegg_edges","kegg_nodes"),file = paste0(save_dir0,"/kegg_info.RData"))
  }
  kegg_info <- list(kegg_edges,kegg_nodes)
}


# translate graph_list to PFPRefnet class
load_PFPRefnet  <- function(graph_list,dict){
  PFPRefnet
}
  
