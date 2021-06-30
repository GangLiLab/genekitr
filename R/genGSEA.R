##' GSEA for a genelist with logFC
##'
##' @param genelist order ranked genelist.
##' @param org species "mm" or "hs".
##' @return a dataframe of gene info.
##' @importFrom stringr str_to_title
##' @export
##' @examples
##' \dontrun{
##'
##' }


genGSEA <- function(genelist,
                    org,
                    ...){

  #--- args ---#
  if (!is.sorted(genelist))
    stop("genelist should be a decreasing sorted vector...")









}
