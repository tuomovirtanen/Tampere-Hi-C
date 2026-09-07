genome_plot <- function(chr, pos, y,
                        point_color = "black",
                        ylab = "Value",
                        segments = NULL,       # <-- new optional argument
                        segment_color = "red", # optional styling
                        line_width = 0.8,
                        genome_build = "hg38") {
  suppressPackageStartupMessages({
    library(dplyr)
    library(GenomeInfoDb)
    library(ggplot2)
  })
  
  # Prepare input
  df <- data.frame(
    chr = chr,
    pos = pos,
    value = y
  )

  # Get chromosome lengths
  chrom_lengths <- seqlengths(Seqinfo(genome = genome_build))
  chrom_lengths <- chrom_lengths[grepl("^chr[0-9XY]+$", names(chrom_lengths))]

  # Chromosome offsets
  chrom_df <- data.frame(
    chr = names(chrom_lengths),
    chr_len = as.numeric(chrom_lengths)
  ) %>%
    mutate(chr_start = c(0, cumsum(chr_len)[-length(chr_len)]),
           chr_end = cumsum(chr_len),
           chr_mid = (chr_start + chr_end) / 2)

  # Merge log2ratio data with genome-wide coordinates
  df_wide <- df %>%
    inner_join(chrom_df, by = "chr") %>%
    mutate(genome_pos = chr_start + pos)

  # Base plot
  gg <- ggplot() +
    geom_point(
      data = df_wide,
      aes(x = genome_pos, y = value),
      color = point_color, size = 0.1
    ) +
    geom_vline(
      data = chrom_df[-nrow(chrom_df), ],
      aes(xintercept = chr_end),
      color = "grey70", linetype = "dashed"
    ) +
    scale_x_continuous(
      name = "Genomic position",
      breaks = chrom_df$chr_mid,
      labels = chrom_df$chr,
      expand = c(0, 0)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    ylab(ylab)

  # --- OPTIONAL: add segment lines ---
  if (!is.null(segments)) {
    seg_wide <- segments %>%
      inner_join(chrom_df, by = c("chr")) %>%
      mutate(
        genome_start = chr_start + start,
        genome_end = chr_start + end
      )
    
    gg <- gg +
      geom_segment(
        data = seg_wide,
        aes(x = genome_start, xend = genome_end,
            y = log2ratio, yend = log2ratio),
        color = segment_color,
        linewidth = line_width
      )
  }

  return(gg)
}

# genome_plot <- function(chr, pos, y,
#                         point_color = "black",
#                         main = "Genome-wide plot",
#                         ylab = "Value") {
#   suppressPackageStartupMessages({
#     library(dplyr)
#     library(GenomeInfoDb)
#     library(ggplot2)
#   })
  
#     df <- data.frame(
#     chr = chr,
#     pos = pos,
#     value = y  # something to plot on y-axis,
#     )

#     # 1. Get chromosome lengths (human hg38 for example)
#     chrom_lengths <- seqlengths(Seqinfo(genome = "hg38"))
#     chrom_lengths <- chrom_lengths[grepl("^chr[0-9XY]+$", names(chrom_lengths))]

# # 2. Make chromosome offset table
#     chrom_df <- data.frame(
#     chr = names(chrom_lengths),
#     chr_len = as.numeric(chrom_lengths)
#     ) %>%
#     mutate(chr_start = c(0, cumsum(chr_len)[-length(chr_len)]),
#             chr_end = cumsum(chr_len),
#             chr_mid = (chr_start + chr_end) / 2)

#     # 3. Join data + compute genome-wide coordinate
#     df_wide <- df %>%
#     inner_join(chrom_df, by = "chr") %>%
#     mutate(genome_pos = chr_start + pos)


# gg <- ggplot() +
#   # First plot FALSE points (bottom layer)
#   geom_point(
#     data = df_wide,
#     aes(x = genome_pos, y = value),
#     color = "black", size = 0.25
#   ) +
#   # Optional vertical chromosome lines
#   geom_vline(
#     data = chrom_df[-nrow(chrom_df), ],
#     aes(xintercept = chr_end),
#     color = "grey70", linetype = "dashed"
#   ) +
#   scale_x_continuous(
#     name = "Genomic position",
#     breaks = chrom_df$chr_mid,
#     labels = chrom_df$chr
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, vjust = 0.5),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank()
#   )+ylab(ylab)


#   return(gg)
# }