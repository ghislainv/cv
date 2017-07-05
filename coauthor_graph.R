if (!"igraph"  %in% installed.packages()[, 1]) {install.packages('igraph')}
if (!"gsubfn"  %in% installed.packages()[, 1]) {install.packages('gsubfn')}
if (!"scholar" %in% installed.packages()[, 1]) {install.packages('scholar')}
require(igraph)
require(gsubfn)
require(scholar)

# make a network graph of Authors
plot_authors <- function(focal_author, prune = 1){
  cat('getting data...\n')
  dat  <- get_publications(focal_author)
  prof <- get_profile(focal_author)
  Auth <- unlist(strsplit(as.vector(dat$author), split = ", "))
  main_auth <- names(sort(table(Auth), T)[1]) # pulls out the most-frequent author. Should always return the name of the focal_author (but formatted in Scholar style)
  complete_list <- sub('...', main_auth, dat$author, fixed = T)
  complete_list <- sub(' jr', '_Jr', complete_list, ignore.case = T)
  complete_list <- complete_list[grepl(sub(".* ", "", main_auth), complete_list)] # this drops out any papers that do NOT have the main author on them, for whatever reason
  Auth <- unlist(strsplit(as.vector(complete_list), split = ", "))
  Auth <- gsubfn(paste(names(unwanted_array),collapse='|'), unwanted_array, Auth) # cleans out all accents
  complete_list <- gsubfn(paste(names(unwanted_array),collapse='|'), unwanted_array, complete_list)
  complete_list <-sapply(complete_list , tolower)
  auth_list <- table(tolower(sub(".* ", "",Auth)))
  
  # for every combination of authors, this indicates on how many papers have they shared authorship
  auth_mat <- matrix(0, nrow = length(auth_list), ncol = length(auth_list))
  colnames(auth_mat) <- names(auth_list)
  for(i in 1:length(auth_list)){
    which.i <- grepl(names(auth_list)[i], complete_list)
    for(j in 1:length(auth_list)){
      auth_mat[i,j] <- sum(which.i & grepl(names(auth_list)[j], complete_list))
    }
  }
  cat('making graph...\n')
  
  #str(auth_mat)
  auth_mat[auth_mat < prune] <- 0 # any co-author who has partipated in fewer than prune papers gets dropped from the graph. 
  ig <- graph_from_adjacency_matrix(auth_mat, mode = 'undirected', weighted = T, diag = F)
  ig <- delete.isolates(ig)
  weights <- E(ig)$weight-prune+1 # downscale the line widths
  COL <- sprintf('gray%1.0f',91-90*weights/max(weights))
  
  # ig
  cat('finding clusters...\n')
  wc <- cluster_spinglass(ig) # this decides on clusters, based on linkage
  
  # and finally plot
  op <- par(mfrow =c(1, 1), las = 1, bty = 'n', lend = 2, mar =c(0, 0,0,0))
  plot.new()
  plot(wc, ig, vertex.label = simpleCap(V(ig)$name), vertex.label.cex = 0.6, vertex.frame.color = NA, edge.width = weights, edge.color = COL, vertex.shape = 'rectangle', vertex.size = pmin(20, strwidth(V(ig)$name)*150), vertex.size2 = 8)
  
  #	plot(ig, mark.border = NA, vertex.color = 'white', vertex.label = simpleCap(V(ig)$name), vertex.label.cex = 0.6, edge.width = weights, edge.color = COL) # use this when defining groups fails becauuse of an un-connected graph
  mtext(prof$name, adj = 0, line = -1, font = 2)
  mtext(prof$affiliation, adj = 0, line = -2)
  legend('topright', lwd = range(weights), title = 'N joint publications', legend = range(E(ig)$weight), seg.len = 4, col = c('gray75', 'gray1'), inset = 0.05)
  par(op)
}


unwanted_array = list('À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='A', 'Å'='A', 'Æ'='A', 'Þ'='B', 'Ç'='C', 'È'='E', 'É'='E', 'Ê'='E', 'Ë'='E', 'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 'Ñ'='N', 'Ò'='O', 'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 'Š'='S', 'š'='s', 'ß'='Ss', 'Ù'='U', 'Ú'='U', 'Û'='U', 'Ü'='U', 'Ý'='Y', 'Ž'='Z', 'ž'='z', 'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='a', 'å'='a', 'æ'='a', 'þ'='b', 'ç'='c', 'è'='e', 'é'='e', 'ê'='e', 'ë'='e', 'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ñ'='n', 'ð'='o', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o', 'ö'='o', 'ø'='o', 'ù'='u', 'ú'='u', 'û'='u', 'ü'= 'u', 'ý'='y', 'ý'='y', 'ÿ'='y')


simpleCap <- function(x) {
  x <- tolower(x)
  paste(toupper(substring(x, 1, 1)), substring(x, 2), sep = "")
}
delete.isolates <- function(graph) {
  isolates <- which(degree(graph) == 0)
  delete.vertices(graph, isolates)
}

gv <- "fv71B7UAAAAJ" # Google Scholar ID
plot_authors(gv)
