#by Karsten Bach 
###This function counts the number of edges between two celltypes in the tissue and normalizes for the number of cells
getContacts <- function(input,rad,knn,expected) {
  nn <- nn2(input[,c("X","Y")],k=knn,radius=rad,searchtype="radius")
  
  indx.knn <- nn[["nn.idx"]]
  m.adj <- Matrix(0, nrow=nrow(indx.knn), ncol=nrow(indx.knn), sparse=TRUE) 
  rownames(m.adj) <- colnames(m.adj) <- rownames(indx.knn) <- rownames(input) <- input$Barcode
  
  for (i in seq_len(nrow(m.adj))) {
    m.adj[i,rownames(indx.knn)[indx.knn[i,]]] <- 1
  }
  diag(m.adj) <- 0 # Remove the self connection
  igr <- graph_from_adjacency_matrix(m.adj)
  igr <- set_vertex_attr(igr,name="annotation",value = unlist(input[names(V(igr)),"annotation"]))
  ct.abundance <- table(input$annotation)
  out <- data.frame()
  for (ct in expected) {
    vs <- V(igr)$name[V(igr)$annotation==ct]
    egds <- unlist(lapply(adjacent_vertices(igr,v=vs),names))
    if (length(egds)==0) { # In case a cell type as 0 edges
      tmp <- data.frame("Edges"=0,
                        "From"=ct,
                        "To"=expected)
    } else {
      tbbl <- table(input[egds,"annotation"])
      normfac <- (length(vs) + ct.abundance[names(tbbl)])
      normcounts <- as.numeric(unname(tbbl)/normfac) * 100
      tmp <- data.frame("Edges"=normcounts,
                        "From"=ct,
                        "To"=names(tbbl))
    }
    mssng <- setdiff(expected,tmp$To)
    if (length(mssng)>0) {
      add  <- data.frame("From"=ct,
                         "To"=mssng,
                         "Edges"=0)
      tmp <- rbind(tmp,add)
    }
    out <- rbind(out,tmp)
  }
  return(out)
}
# Get the difference in edges between two
getDiff <- function(input,diffFrom,diffTo) {
  input$Area <- factor(input$Area,levels=c(diffFrom,diffTo))
  tmp <- split(input,f=input$Area)
  tmp <- lapply(tmp, function(DF) {
    rownames(DF) <- DF$Label <- paste0(DF$From,".",DF$To,".",DF$Slide)
    return(DF)
  })
  tmp[[2]] <- tmp[[2]][rownames(tmp[[1]]),]
  out <- tmp[[1]]
  out[,paste0("EdgeDiff_To_",diffTo)] <- tmp[[1]][,"Edges"] - tmp[[2]][,"Edges"]
  return(out)
}


runTest <- function(input, con, k, rad, n_it, BPPARAM, expected, diffFrom, diffTo) {
  #Run shuffeling
  chance <- bplapply(1:n_it, BPPARAM=BPPARAM, function(i) {
    # 1st Shuffle CellType assignment
    input.rnd <- input
    input.rnd$Barcode <- sample(input.rnd$Barcode, length(input.rnd$Barcode), replace=FALSE)
    input.rnd$annotation <- input[input.rnd$Barcode,"annotation"]
    input.rnd  <- split(input.rnd, f=input$Area)
    # 2nd Compute statistics
    boot.list <- lapply(input.rnd,function(INPUT) {
      cons <- getContacts(INPUT,k=k,rad=rad,expected=expected)
      cons$Area <- unique(INPUT$Area)
      cons$Slide <- unique(INPUT$Slide)
      return(cons)
    })
    boot.con <- do.call(rbind,boot.list)
    boot.con <- getDiff(input=boot.con,diffFrom=diffFrom,diffTo=diffTo)
    boot.con$Iteration <- paste0("Iteration_",i)
    return(data.frame(boot.con))
  })
  
  chance <- do.call(rbind,chance)
  
  con <- data.frame(con) # Protect yourself from stupid tibble
  colname <- paste0("EdgeDiff_To_",diffTo)
  out <- bplapply(unique(con$Label), BPPARAM=BPPARAM, FUN=function(LABEL) {
    measured <- con[con$Label==LABEL,colname]
    emp.dist <- chance[chance$Label==LABEL,colname]
    pval <- min(c(sum(measured >= emp.dist) + 1,sum(measured <= emp.dist) + 1) / (length(emp.dist) + 1))
    Zscore <- (measured - mean(emp.dist)) / sd(emp.dist)
    out <- con[con$Label==LABEL,]
    out$Pval[out$Label==LABEL] <- pval
    out$Zscore[out$Label==LABEL] <- Zscore
    #		 out$Dist[out$Label==LABEL] <- list(emp.dist)
    return(data.frame(out))
  })
  out <- do.call(rbind,out)
}

cellCellContactMap <- function(pmat, order = NULL, exclude = NULL) {
  # this function takes output from cellCellContact() 
  # and gives a ggplot object of the graph
  # out = mat_p_sym
  # if order not given then perform hclust to get ordering
  require(reshape)
  require(ggplot2)
  mat_p_df = melt(pmat)
  colnames(mat_p_df) <- c("CellType_1","CellType_2", "NInteractions")
  if (is.null(order)) {
    hc = hclust(dist(pmat))
    mat_p_df$CellType_1 <- factor(mat_p_df$CellType_1, levels = 
                                    hc$labels[hc$order])
    mat_p_df$CellType_2 <- factor(mat_p_df$CellType_2, levels = 
                                    hc$labels[hc$order])
  } else {
    mat_p_df$CellType_1 <- factor(mat_p_df$CellType_1, levels = 
                                    order)
    mat_p_df$CellType_2 <- factor(mat_p_df$CellType_2, levels = 
                                    order)
  }
  mat_p_df$keep = as.numeric(mat_p_df$CellType_1) >= as.numeric(mat_p_df$CellType_2) 
  g = ggplot(subset(mat_p_df, keep & 
                      (!CellType_1 %in% exclude) & 
                      (!CellType_2 %in% exclude)),
             aes(x = CellType_1, y = CellType_2, fill =NInteractions)) + 
    geom_tile() +
    # geom_text(aes(label = Sig, y = CellType_2 - 0.025), size = 20) +
    # geom_point(pch = "*", size = 10, data = subset(mat_p_df, Sig == "*" & keep)) +
    theme_classic() +
    theme(axis.line = element_blank()) +
    theme(axis.ticks = element_blank()) +
    # scale_y_continuous(position = "right") +
    scale_y_discrete(position = "right") +
    theme(axis.text = element_text(size = 14)) +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_fixed() +
    # guides(fill = guide_legend(title = "")) +
    scale_fill_gradient2(high ="#9E0142",
                         mid = "white",
                         low = "#5182BB",
                         # na.value = "white",
                         midpoint = 0, #limits = c(0,1),
                         limits = c(-1,1),
                         # labels = c("Segregated", "Integrated")
    ) +
    theme(legend.position = "top") + 
    theme(legend.text = element_text(size = 15)) +
    theme(legend.key.width = unit(1, "in")) +
    guides(fill = guide_colourbar(title.position = "top",
                                  title = "",
                                  title.hjust = 0.5,
                                  ticks = FALSE,
                                  reverse = TRUE)) +
    NULL
  # print(g)
  return(g)
}
