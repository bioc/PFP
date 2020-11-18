# kegg_download
#' @title download kegg KGML files
#' @description This function will download all kegg KGML files assigned by
#'  \code{spec}.
#' @param spec, a character, refers to the species names in kegg, such as "hsa",
#' "mmu"...
#' @param file_root, a character,refers to the root you want to save kegg
#' pathway kgml files in.
#' @details Downloading all kegg KGML files assigned by \code{spec} from
#'  https://www.kegg.jp/kegg/xml/,
#'  which may take tens of minutes.
#' @return the kegg KGML files
#' @examples
#' # download the kegg network
#' kegg_download(spec = "hsa", file_root="~/Desktop")
#' @export
kegg_download <- function(spec,file_root="."){
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
  download.file(url = url_pathway_names,
                destfile = paste0(file_root,"/",spec,"_pathways.txt"))

  file1 = read.csv(file = paste0(file_root,"/",spec,"_pathways.txt"),
                   header = FALSE,
                   sep = "\t")

  if (length(unique(substr(as.vector(file1$V1),1,5)))!=1){
    stop("You downloaded the wrong pathway info file! Please check your network
         connection!")
  }
  name <- substr(as.vector(file1$V1),6,13)
  # repeat several times for complete download
  download_repeat <- function(name,download_dir){
    loaded_files <- list.files(download_dir)
    name_undownload <- name[!vapply(X = name,
                                    FUN = function(x)(
                                      paste0(x,".xml") %in% loaded_files),
                                    TRUE)]
    error_files <- vapply(name_undownload,
                          function(x)(
                            "try-error" %in% class(
                              download.file(
                                url = paste0("http://rest.kegg.jp/get/",x,
                                             "/kgml"),
                                destfile = paste0(download_dir,"/",x,".xml"),
                                quiet=TRUE))),
                          TRUE)
    if (length(name_undownload)!=0){
      name_undownload <- name_undownload[error_files]
    }
    name_undownload
  }

  print(paste0("Downloading ",spec," KEGG pathays..."))
  print("This may take you ten minutes, depending on your network.")
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
    print(paste0("You have completely downloaded all ",spec,
                 " kgml files in kegg pathways!"))
  }
}



# trans_xml2graph
#' @title translate kgml files to graphNEl
#' @description This function will translate all kegg KGML files in path
#' \code{file_dir}.
#' @param file_dir, a character, refers to the file_path where kegg KGML files
#' are stored.
#' @details translate all kegg KGML files in path \code{file_dir}. It will
#' return a list of \code{graphNEL}
#' @return Trans the xml to graph
#' @examples
#' # Download the kegg
#' #kegg_download(spec = "hsa", file_root="~/Desktop")
#' # Trans the xml to graph
#' #graph_list <- trans_xml2graph("~/Desktop")
#' # Show the graph
#' data(gene_list_hsa)
#' gene_list_hsa
#' @export
trans_xml2graph <- function(file_dir){
  if (substr(file_dir,nchar(file_dir),nchar(file_dir)) == "/"){
    file_dir <- substr(file_dir,1,(nchar(file_dir)-1))
  }
  xml_list <- list.files(file_dir)
  xml_select <- unlist(data.frame(
    strsplit(x=xml_list,split = "\\."))[2,]) == "xml"
  xml_list <- xml_list[xml_select]
  if (length(xml_list)==0){
    stop("Please input a file dir which contains at least one kgml file!")
  }
  graph_list <- vapply(X = xml_list, FUN = function(file_name) try(
    parseKGML2Graph(paste0(file_dir,"/",file_name))))

  # remove pathways which failed to be translated
  right_trans <- vapply(graph_list,function(x)is(x,"graphNEL"),TRUE)
  wrong_trans_file <- names(graph_list)[!right_trans]
  if (length(wrong_trans_file)!=0){
    print("These files failed to translate! You can manually redownload them
          in KEGG website.")
    print(wrong_trans_file)
  }
  graph_list <- graph_list[right_trans]
  # remove pathways which have no elements in nodes or edges
  zero_id <- vapply(X = graph_list,FUN = function(x)(
    length(nodes(x))==0 & length(edgeData(x))==0),
    TRUE)
  graph_list <- graph_list[!zero_id]
  names(graph_list) <- t(data.frame(strsplit(x = names(graph_list),
                                             split = "\\.")))[,1]

  spec <- substr(names(graph_list)[1],1,3)
  # remove the specices label in ENTREZID  gene id
  fun_trans <- function(id,graph_list){
    graph0 <- graph_list[[id]]
    nodes(graph0) <- as.vector(vapply(nodes(graph0),
                                      function(x)gsub(
                                        pattern = paste0(spec,":"),
                                        replacement = "",
                                        x = x),
                                      character(1)))
    graph1 <- new(Class = "graphNEL",
                  nodes=nodes(graph0),
                  edgeL=edgeL(graph0),
                  edgemode='directed')
  }
  graph_list <- vapply(names(graph_list), FUN = fun_trans,graph_list)
  # save(list = c("graph_list"),file = paste0(file_dir,"/graph_list.RData"))
  # graph_list
}




# trans_xml2graph
#' @title translate graph_list to PFPRefnet class
#' @description This function will translate all graphs in \code{graph_list} to
#' a \code{\link{PFPRefnet-class}} object.
#' @param graph_list, a list of \code{\link{graphNEL}}.
#' @param pathway_info, a data.frame, which contains all kegg pathways "index",
#' "id","name","group","species"
#' @details translating all graphs in \code{graph_list} to a
#' \code{\link{PFPRefnet-class}} object.
#' The pathway_info can be designed by yourself, but the colnames must be
#' "index","id","name","group" and "species".
#' @return a PFPRefnet
#' @examples
#' # Load the info of the pathway
#' data(pathway_info)
#' pathway_info
#' @export
trans_graph2PFPRefnet  <- function(graph_list,pathway_info){
  spec <- substr(names(graph_list)[1],1,3)
  name_len <- nchar(names(graph_list)[1])
  id_graphlist <- data.frame(id = substr(names(graph_list),4,name_len))
  pathway_info <- merge(id_graphlist,pathway_info,by="id",all.x=TRUE)
  pathway_info[["id"]] <- paste0(rep(spec,nrow(pathway_info)),pathway_info$id)
  pathway_info[["species"]] <- rep(spec,nrow(pathway_info))
  pathway_info <- pathway_info[order(pathway_info$index,decreasing = FALSE),
                               c("index","id","name","group","species")]
  graph_list <- graph_list[pathway_info$id]
  PFPRefnet <- new(Class = "PFPRefnet",
                   network = graph_list,
                   net_info = pathway_info)
}


# library(KEGGgraph) # parseKGML2Graph
# library(graph) # nodes edgeL edgeData
# #library(igraph)
# #library(BioNet)
# test
# spec <- "mmu"
# kegg_download(spec,file_root="/home/zx/文档/test")
# file_dir <- "/home/zx/文档/test/kgml/mmu/"
# graph_list <- trans_xml2graph(file_dir)
# load("/home/zx/文档/PFP/RData/pathway_info.RData")
# PFPRefnet_mmu <- trans_graph2PFPRefnet(graph_list,pathway_info)
# save(list = c("PFPRefnet_mmu"),file = "/home/zx/文档/kangqichuang/PFP/data/PFPRefnet_mmu.RData")
