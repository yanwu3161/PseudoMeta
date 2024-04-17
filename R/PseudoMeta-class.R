#' @export PseudoMeta
setClass("PseudoMeta", slots = c(
  assays = "list",
  analysis = "list",
  reductions = "list",
  meta.data = "data.frame",
  ID = "character"
))