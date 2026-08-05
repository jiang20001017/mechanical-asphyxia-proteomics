library(ggplot2)
library(ggrepel)

differential_results <- read.csv(
    "output/dma_hs_statistics.csv",
    check.names = FALSE
)

differential_results$change <- ifelse(
    differential_results$Sig == "Up",
    "Up",
    ifelse(differential_results$Sig == "Down", "Down", "Stable")
)
differential_results$change <- factor(
    differential_results$change,
    levels = c("Down", "Stable", "Up")
)
print(table(differential_results$change))

labelled_genes <- differential_results[
    differential_results$Gene_symbol %in% c("NDUFS8", "SUCLG1"),
]

volcano_plot <- ggplot(
    differential_results,
    aes(x = logFC, y = -log10(pvalue), color = change)
) +
    geom_point(alpha = 0.7, size = 0.8) +
    scale_color_manual(
        values = c(Down = "#005696", Stable = "grey", Up = "#AA0000")
    ) +
    geom_vline(
        xintercept = c(-0.263, 0.263),
        lty = 2,
        col = "grey40",
        lwd = 0.4
    ) +
    geom_hline(
        yintercept = -log10(0.05),
        lty = 2,
        col = "grey40",
        lwd = 0.4
    ) +
    labs(x = "Log2(Fold Change)", y = "-Log10(P.Value)") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        aspect.ratio = 1
    ) +
    scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, 2)) +
    scale_y_continuous(limits = c(-0.5, 5)) +
    ggrepel::geom_text_repel(
        data = labelled_genes,
        aes(label = Gene_symbol),
        color = "black",
        size = 3,
        box.padding = 1,
        point.padding = 0.5,
        segment.color = "black",
        segment.size = 0.3,
        min.segment.length = 0,
        force = 6,
        max.overlaps = Inf,
        show.legend = FALSE
    )

volcano_plot

ggsave("output/dma_hs_volcano_plot.pdf", volcano_plot, width = 5, height = 5)
