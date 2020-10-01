# id translate
#' @title translate other id type to ENTREZID
#' @description translate other id type to ENTREZID
#' @param gene_list, a vector of characters
#' @param gene_id_type, a character,the type of gene ID, "ENSEMBL","GO","SYMBOL" and so on.
#' @param drop, a logical,whether to drop NA when translating id type
#' @details, translate other id type to ENTREZID,
#' @examples
#' \dontrun{
#' library(org.Mm.eg.db)
#' data(gene_list)
#' gene_n <- trans_id(gene_list,gene_id_type)
#' }
#' @export
trans_id <- function(gene_list,gene_id_type,gene_info_db,drop=FALSE){
  gene_t <- gene_list
  if (gene_id_type != "ENTREZID"){
    gene_n <- bitr(gene_t,fromType = gene_id_type,
                   toType = "ENTREZID",#"ENSEMBL"
                   OrgDb = gene_info_db,drop = drop)
    gene_n <- gene_n[!duplicated(gene_n$ENSEMBL),]
  }else{
    gene_n <- data.frame(ENTREZID=gene_t)
  }
  return (gene_n)
}
