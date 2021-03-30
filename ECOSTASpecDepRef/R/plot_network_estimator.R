normalize_weight <- function(metric) {
  metric0 <- metric - diag(diag(metric))
  weight <- (metric0 - min(metric0)) / (max(metric0) - min(metric0))
  weight
}

adjust_df_names <- function(anames) {
  #  Expected input
  # anames <- c("end", "start", "weight", "band", "band_id", "metric", "xpos.x", "ypos.x", "xpos.y", "ypos.y")
  #  Expected output
  # names(network_data) = c("end", "start", "weight", "band", "band_id", "metric", "x", "y", "xend", "yend")
  anames[anames == "xpos.x"] <- "x"
  anames[anames == "ypos.x"] <- "y"
  anames[anames == "xpos.y"] <- "xend"
  anames[anames == "ypos.y"] <- "yend"
  anames
}

stretch_edge_distances <- function(network_data, arrow_distance_x =0.035, arrow_distance_y=0.035) {
  dx <- (network_data$xend - network_data$x) * arrow_distance_x
  network_data$x <- network_data$x + dx
  network_data$xend <- network_data$xend - dx
  dy <- (network_data$yend - network_data$y) * arrow_distance_y
  network_data$y <- network_data$y + dy
  network_data$yend <- network_data$yend - dy
  network_data
}

create_network_layout <- function(nodes, radius=0.5, xcenter=0.5, ycenter=0.5) {
  node_info <- list(xpos=c(), ypos=c(), node=nodes, node_id=nodes)
  angle <- 360 / length(nodes)
  angle0 <- 360 / length(nodes) * -0.5
  for(i in 1:length(nodes)) {
    x <- xcenter + radius * cos((angle0 + angle * i) * pi / 180)
    y <- xcenter + radius * sin((angle0 + angle * i) * pi / 180)
    node_info$xpos[i] <- x
    node_info$ypos[i] <- y
  }
  as.data.frame(node_info)
}

#' @export
plot_network_estimator <- function(X, fs, network_estimator, is_symmetric=FALSE, label_id=1, selected_bands=c("delta", "alpha", "gamma"), nodes=NULL, colors=c("steelblue", "darkorange", "olivedrab3", "firebrick3"), padding=0.15, node_radius=10, color_edge="grey50", threshold_rel_max=0.1, metric_name="|PDC|", return_metric=FALSE){
  #
  #
  if(as.numeric(label_id)) {
    label_id <- intToUtf8(utf8ToInt("A") + label_id - 1)
  }
  if(is.null(nodes)) {
    nodes <- paste0("X", 1:ncol(X))
  }
  stopifnot(length(colors) >= ncol(X))
  stopifnot(length(colors) >= length(nodes))
  N <- nrow(X)
  names(colors) <- nodes
  node_info <- create_network_layout(nodes)
  x_min <- min(node_info$xpos)
  x_max <- max(node_info$xpos)
  y_min <- min(node_info$ypos)
  y_max <- max(node_info$ypos)
  #
  edge_info <- list(start=c(), end=c(), weight=c(), band=c(), band_id=c())
  m <- 1
  for(k in 1:length(selected_bands)){
    metric <- network_estimator(
      X,
      fs=fs,
      selected_band=selected_bands[[k]]
    )
    weight <- metric
    #weight <- normalize_weight(metric)
    for(i in 1:ncol(X)){
      for(j in 1:ncol(X)){
        if(i == j) next;
        if(is_symmetric && i <= j) next; #Coherence is symmetric
        edge_info$band[m] <- selected_bands[k]
        edge_info$band_id[m] <- (k) * (as.numeric(i > j) * 2 - 1)
        edge_info$start[m] <- nodes[i]
        edge_info$end[m] <- nodes[j]
        edge_info$weight[m] <- weight[i, j]
        edge_info$metric[m] <- metric[i, j]
        m <- m + 1
      }
    }
  }
  
  edge_info <- as.data.frame(edge_info)
  edge_info$weight <- abs(edge_info$weight)
  network_data <- merge(
    merge(edge_info, node_info, by.x="start", by.y="node_id"),
    node_info, by.x="end", by.y="node_id")
  network_data$node.x = NULL
  network_data$node.y = NULL
  names(network_data) <- adjust_df_names(names(network_data))
  network_data$band <- as.character(network_data$band)
  network_data$curvature <- (network_data$band_id)
  network_data <- stretch_edge_distances(network_data)
  threshold <- threshold_rel_max * max(network_data$weight) # min(v) + 0.4 *( max(v) - min(v) )
  threshold <- threshold_rel_max
  print(sprintf("Threshold: %.5f", threshold))
  print(network_data)
  raw_network_data <- network_data
  network_data <- network_data[network_data$weight > threshold,]
  
  p <- ggplot(data=network_data) +
    geom_point(data=node_info, aes(x=xpos, y=ypos, colour=node), size=node_radius) +
      lapply(network_data$curvature, function(curvature){
        geom_curve(
          aes(
            x=x, y=y, xend=xend, yend=yend,
            size=weight, linetype=band
          ),
          #lineend="round",
          #linejoin="round",
          arrow=ggplot2::arrow(angle=40, length=ggplot2::unit(if(is_symmetric) 0 else 0.15, "inches"), ends="last", type="open"), 
          color=if(is_symmetric) (color_edge) else (if(curvature >= 0) "darkred" else "darkblue"),
          angle=90 - 5 * curvature,
          curvature=abs(curvature) * 0.1 + 0,
          data=network_data[network_data$curvature==curvature,]
        )
      }) +
      scale_size(name=metric_name, range=c(0.1, 2), limits=c(threshold, 1), trans="exp") + 
      #scale_size(name=metric_name, range=c(0.1, 2), limits=c(threshold, 1), trans="identity") + 
      geom_point(data=node_info, aes(x=xpos, y=ypos, colour=node), size=node_radius, alpha=0.7) +
      geom_text(data=node_info, aes(x=xpos, y=ypos, label=node_id), color="white") +
      scale_colour_manual(values=colors) +
      scale_linetype_manual(values=c("solid", "longdash", "twodash", "dotted", "dashed")) +
      coord_cartesian(xlim=c(
        x_min - padding * (x_max - x_min),
        x_max + padding * (x_max - x_min)
      ), ylim=c(
        y_min - padding * (y_max - y_min),
        y_max + padding * (y_max - y_min)
      )) +
      labs(title=sprintf("%s) Network links", label_id)) +
      csda_theme(plain=TRUE) +
      theme(legend.position="left")
    if(return_metric) list(plot=p, data=raw_network_data)
    else p
}
