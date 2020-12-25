# kegg_download
#' @title download kegg KGML files
#' @description This function will download all kegg KGML files assigned by
#'  \code{spec}.
#' @param spec, a character, refers to the species names in kegg, such as "hsa",
#' "mmu"...
#' @param file_root, a character,refers to the root you want to save kegg
#' pathway kgml files in.
#' @param test_mode, please set whether to test this function.
#' @details Downloading all kegg KGML files assigned by \code{spec} from
#'  https://www.kegg.jp/kegg/xml/,
#'  which may take tens of minutes.
#' @return the kegg KGML files
#' @examples
#' kegg_download(spec,file_root=".", test_mode=TRUE)
#' @export
kegg_download <- function(spec,file_root=".",test_mode=FALSE){
  if(test_mode){
    print("You have started test mode, please change test mode to FALSE ")
  }
  else{
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
                                  quiet=TRUE))),TRUE)
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

}



# trans_xml2graph
#' @title translate kgml files to graphNEl
#' @description This function will translate all kegg KGML files in path
#' \code{file_dir}.
#' @param file_dir, a character, refers to the file_path where kegg KGML files
#' are stored.
#' @param test_mode, please set whether to test this function.
#' @details transform all KEGG KGML files downloaded by the function 
#' kegg_download() in path \code{file_dir} to the graphNEL object
#' @return a list of \code{graphNEL}
#' @examples
#' trans_xml2graph(file_dir, test=TRUE)
#' @export
trans_xml2graph <- function(file_dir,test_mode=FALSE){
  if(test_mode){
    print("You have started test mode, please change test mode to FALSE ")
  }
  else{
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
    graph_list <- lapply(X = xml_list, FUN = function(file_name) try(
      parseKGML2Graph(paste0(file_dir,"/",file_name))))
    names(graph_list) <- xml_list
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
    graph_list2 <- lapply(names(graph_list), FUN = fun_trans,graph_list)
    names(graph_list2) <- names(graph_list)
    graph_list2
    # save(list = c("graph_list"),file = paste0(file_dir,"/graph_list.RData"))
    # graph_list
  }
}

# get_pathway_info
#' @title get pathway info of a species in KEGG
#' @description This function helps get pathway info of a species in KEGG.
#' @param spec, a character, refers to the species in KEGG. hsa, mmu...
#' @details, get pathway info of a species in KEGG. It will
#' return a data.frame.
#' @return a data.frame whose colnames contains "index","id","name" and "group"
#' @examples
#' pathway_info <- get_pathway_info("hsa")
#' @export
get_pathway_info <- function(spec){
  url <- paste0("https://www.kegg.jp/kegg-bin/show_organism",
                "?menu_type=pathway_maps&org=",spec)
  page <- readLines(url)
  group_num <- grep("^<b>.*</b>$", page)
  pathway_num <- grep("^0.*</a><br>$", page)
  group_num2 <- c(group_num,length(page))
  
  list_pathway_num <- lapply(seq(length(group_num)),
                             function(i)intersect(group_num2[i]:group_num2[i+1],
                                                  pathway_num))
  
  pathway_info <- lapply(seq(length(group_num)),function(i){
    group0 <- substr(page[group_num[i]],4,nchar(page[group_num[i]])-4)
    pathway_info0 <- lapply(X = list_pathway_num[[i]],
                            FUN = function(j){
                              str_sub <- sub(pattern = "[0-9].*pathway\\?",
                                             replacement = "",x = page[j])
                              c(substr(str_sub,1,8),
                                substr(str_sub,11,nchar(str_sub)-8),
                                group0)
                            })
    pathway_info0 <- data.frame(t(data.frame(pathway_info0)))
    colnames(pathway_info0) <- c("id","name","group")
    rownames(pathway_info0) <- seq_len(1):nrow(pathway_info0)
    pathway_info0
    })
  pathway_info <- do.call(rbind,pathway_info)
  pathway_info <- cbind(data.frame(index=seq_len(1):nrow(pathway_info)),pathway_info)
}




# trans_graph2PFPRefnet
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
#' data(PFPRefnet_hsa)
#' PFPRefnet_hsa
#' @export
trans_graph2PFPRefnet  <- function(graph_list,pathway_info){
  spec <- substr(names(graph_list)[1],1,3)
  name_len <- nchar(names(graph_list)[1])
  id_graphlist <- data.frame(id = names(graph_list))
  pathway_info <- merge(id_graphlist,pathway_info,by="id",all.x=TRUE)
  pathway_info[["id"]] <- pathway_info$id
  pathway_info[["species"]] <- rep(spec,nrow(pathway_info))
  pathway_info <- pathway_info[order(pathway_info$index,decreasing = FALSE),
                               c("index","id","name","group","species")]
  graph_list <- graph_list[pathway_info$id]
  PFPRefnet <- new(Class = "PFPRefnet",
                   network = graph_list,
                   net_info = pathway_info)
}


# get_PFPRefnet
#' @title get a PFPRefnet for a species
#' @description This function helps update the latest PFPRefnet odject for a species
#' @param spec, a character, refers to the species in KEGG. hsa, mmu...
#' @param file_root, a character, file dir to download the kgml files.
#' @param test_mode, please set whether to test this function.
#' @details, gupdate the latest PFPRefnet odject for a species in KEGG. It will
#' return a PFPRefnet object.
#' @return a PFPRefnet object.
#' @examples
#' PFPRefnet1 <- get_PFPRefnet("hsa",".",test_mode=TRUE)
#' @export
get_PFPRefnet <- function(spec,file_root=".",test_mode=FALSE){
  if(test_mode){
    print("You have started test mode, please change test mode to FALSE ")
  }else{
    if (substr(file_root,nchar(file_root),nchar(file_root)) == "/"){
      file_root <- substr(file_root,1,(nchar(file_root)-1))
    }
    kegg_download(spec=spec,file_root=file_root)
    graph_list <- trans_xml2graph(file_dir=paste0(file_root,"/kgml/",spec))
    pathway_info <- get_pathway_info(spec = spec)
    trans_graph2PFPRefnet(graph_list=graph_list,pathway_info=pathway_info)
  }
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
