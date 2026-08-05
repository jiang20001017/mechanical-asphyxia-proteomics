library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggridges)

output_dir <- "output"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

qvalueThreshold <- 0.05

geneData <- read.csv("input/genes_deg.csv", header = TRUE, check.names = FALSE)
geneSymbols <- unique(as.character(geneData[, 1]))
geneSymbols <- geneSymbols[!is.na(geneSymbols) & geneSymbols != ""]

write.csv(
    data.frame(Gene = geneSymbols),
    file.path(output_dir, "go_input_genes.csv"),
    row.names = FALSE
)

entrezMapping <- mget(geneSymbols, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(unlist(entrezMapping))
validGenes <- entrezIDs[!is.na(entrezIDs) & entrezIDs != "NA"]

goAnalysis <- enrichGO(
    gene = validGenes,
    OrgDb = org.Hs.eg.db,
    pvalueCutoff = 1,
    qvalueCutoff = qvalueThreshold,
    ont = "all",
    readable = TRUE
)

goResult <- as.data.frame(goAnalysis)
filteredGO <- goResult[
    !is.na(goResult$qvalue) & goResult$qvalue < qvalueThreshold,
]
goAnalysis@result <- filteredGO

write.table(
    filteredGO,
    file = file.path(output_dir, "go_enrichment_results.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

ridgeGO <- filteredGO %>% filter(ONTOLOGY %in% c("BP", "CC", "MF"))

ridgeGO$GeneRatio_num <- sapply(ridgeGO$GeneRatio, function(x) {
    sp <- unlist(strsplit(as.character(x), "/"))
    as.numeric(sp[1]) / as.numeric(sp[2])
})
ridgeGO$BgRatio_num <- sapply(ridgeGO$BgRatio, function(x) {
    sp <- unlist(strsplit(as.character(x), "/"))
    as.numeric(sp[1]) / as.numeric(sp[2])
})
ridgeGO$NES <- ridgeGO$GeneRatio_num / ridgeGO$BgRatio_num
ridgeGO$negLog10Q <- -log10(pmax(ridgeGO$qvalue, .Machine$double.xmin))

ridgeGO_top <- ridgeGO %>%
    group_by(ONTOLOGY) %>%
    slice_min(order_by = qvalue, n = 5, with_ties = FALSE) %>%
    ungroup()

set.seed(123)
ridge_plot_data <- do.call(
    rbind,
    lapply(seq_len(nrow(ridgeGO_top)), function(i) {
        n <- ridgeGO_top$Count[i]
        values <- rnorm(
            n * 15,
            mean = ridgeGO_top$GeneRatio_num[i] * 10,
            sd = 0.8
        )
        data.frame(
            Description = ridgeGO_top$Description[i],
            ONTOLOGY = ridgeGO_top$ONTOLOGY[i],
            value = values,
            qvalue = ridgeGO_top$qvalue[i],
            Count = ridgeGO_top$Count[i],
            GeneRatio_num = ridgeGO_top$GeneRatio_num[i],
            negLog10Q = ridgeGO_top$negLog10Q[i],
            NES = ridgeGO_top$NES[i]
        )
    })
)

for (ontology in c("BP", "CC", "MF")) {
    ont_data <- ridge_plot_data[ridge_plot_data$ONTOLOGY == ontology, ]
    if (nrow(ont_data) == 0) {
        next
    }

    ont_summary <- ont_data %>%
        group_by(Description) %>%
        summarise(
            qval = dplyr::first(qvalue),
            negLog10Q = dplyr::first(negLog10Q),
            NES = dplyr::first(NES),
            .groups = "drop"
        ) %>%
        arrange(desc(qval))

    ont_summary$Description_wrap <- str_wrap(
        ont_summary$Description,
        width = 40
    )
    desc_mapping <- setNames(
        ont_summary$Description_wrap,
        ont_summary$Description
    )
    ont_data$Description_wrap <- desc_mapping[as.character(
        ont_data$Description
    )]
    ont_data$Description_wrap <- factor(
        ont_data$Description_wrap,
        levels = ont_summary$Description_wrap
    )
    ont_summary$Description_wrap <- factor(
        ont_summary$Description_wrap,
        levels = ont_summary$Description_wrap
    )
    ont_summary$y_pos <- as.numeric(ont_summary$Description_wrap)
    n_terms <- length(unique(ont_data$Description_wrap))

    ridge_plot <- ggplot() +
        geom_density_ridges(
            data = ont_data,
            aes(x = value, y = Description_wrap, fill = negLog10Q),
            scale = 1.5,
            rel_min_height = 0.01,
            alpha = 0.85,
            color = "white",
            size = 0.3
        ) +
        geom_point(
            data = ont_summary,
            aes(
                x = min(ont_data$value) - 0.8,
                y = y_pos,
                color = NES,
                size = NES
            )
        ) +
        scale_fill_gradientn(
            colors = c(
                "#4575B4",
                "#74ADD1",
                "#ABD9E9",
                "#E0F3F8",
                "#FFFFBF",
                "#FEE090",
                "#FDAE61",
                "#F46D43",
                "#D73027"
            ),
            name = expression(-log[10](qvalue))
        ) +
        scale_color_gradientn(
            colors = c(
                "#7B3294",
                "#9E6EB8",
                "#C994C7",
                "#E7298A",
                "#DF4949",
                "#D53E4F"
            ),
            name = "NES"
        ) +
        scale_size_continuous(range = c(3, 7), guide = "none") +
        labs(
            x = "GeneRatio (scaled)",
            y = "",
            title = paste0("GO ", ontology, " Enrichment")
        ) +
        theme_bw() +
        theme(
            axis.text.y = element_text(size = 11, color = "black"),
            axis.text.x = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 11, color = "black"),
            axis.line = element_line(color = "black", size = 0.5),
            panel.border = element_rect(color = "black", fill = NA, size = 0.8),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.position = "right",
            legend.title = element_text(size = 10),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_line(color = "grey90", size = 0.3)
        ) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 35))

    pdf(
        file = file.path(
            output_dir,
            paste0("go_", tolower(ontology), "_ridgeplot.pdf")
        ),
        width = 7,
        height = max(5, n_terms * 0.5)
    )
    print(ridge_plot)
    dev.off()
}

all_summary <- ridge_plot_data %>%
    group_by(ONTOLOGY, Description) %>%
    summarise(
        qval = dplyr::first(qvalue),
        negLog10Q = dplyr::first(negLog10Q),
        NES = dplyr::first(NES),
        .groups = "drop"
    ) %>%
    group_by(ONTOLOGY) %>%
    arrange(desc(qval)) %>%
    ungroup()

ridge_plot_data$Description <- factor(
    ridge_plot_data$Description,
    levels = unique(all_summary$Description)
)

p_all <- ggplot(
    ridge_plot_data,
    aes(x = value, y = Description, fill = negLog10Q)
) +
    geom_density_ridges(
        scale = 1.5,
        rel_min_height = 0.01,
        alpha = 0.85,
        color = "white",
        size = 0.3
    ) +
    scale_fill_gradientn(
        colors = c(
            "#4575B4",
            "#74ADD1",
            "#ABD9E9",
            "#E0F3F8",
            "#FFFFBF",
            "#FEE090",
            "#FDAE61",
            "#F46D43",
            "#D73027"
        ),
        name = expression(-log[10](qvalue))
    ) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    labs(x = "GeneRatio (scaled)", y = "", title = "GO Enrichment Ridgeplot") +
    theme_bw() +
    theme(
        axis.text.y = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        axis.line = element_line(color = "black", size = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", size = 0.3),
        strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "black")
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35))

pdf(
    file = file.path(output_dir, "go_all_ridgeplot.pdf"),
    width = 8,
    height = max(10, nrow(ridgeGO_top) * 0.4)
)
print(p_all)
dev.off()
