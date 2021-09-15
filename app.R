# OBmap Shiny App
# Kevin Zhu
# Matsunami Lab, Duke University

# define UI for application that plots 3D voxels with OR glomeruli probabilities

#v1, 2021-09-15
# allows user to
# 1. select an OR of interest and number of top probability voxels to view interactively

#plan for v2 is to add functionality for users to:
# a. view additional OR information such as aliases, class, OE zone, etc.
# b. view predictions for OR of interest
# c. view multiple ORs at once (preset groups will include class, OE zone, TAARs, FIenriched, etc.)
# d. manually draw line of mirror symmetry as well as other constraints and re-predict

#load packages
library(shiny)
library(dqshiny)
library(plotly)
library(tidyverse)

#load model output and get gene names
probabilities <- read_csv("ob3d_posteriors_210808v1.csv") #71MB, too large for github
gene_names <- probabilities %>%
    select(olfrname) %>%
    arrange(olfrname) %>%
    unique() %>%
    pull(olfrname)

#given a gene and a cutoff number of top ranked voxels,
#rank and filter probabilities then plot interactively in 3D
RankPosterior3D <- function(gene, rankx = 2011) {
    plot_p50 <- probabilities %>%
        filter(olfrname == gene) %>%
        mutate(p50_rank = min_rank(desc(p50))) %>%
        filter(p50_rank <= rankx) %>%
        mutate(text = paste('Rank:', p50_rank,
                            '<br>p50:', signif(p50, digits = 3),
                            '<br>AP:', AntPos,
                            '<br>VD:', VenDor,
                            '<br>ML:', MedLat))
    blank_p50 <- probabilities %>%
        filter(olfrname == gene) %>%
        mutate(p50_rank = min_rank(desc(p50))) %>%
        filter(p50_rank > rankx) %>%
        mutate(text = paste('Rank:', p50_rank,
                            '<br>p50:', signif(p50, digits = 3),
                            '<br>AP:', AntPos,
                            '<br>VD:', VenDor,
                            '<br>ML:', MedLat))

    rank_plot <- plot_ly(type = "scatter3d", mode = "markers") %>%
        add_trace(data = blank_p50,
                  x = ~AntPos, y = ~MedLat, z = ~VenDor,
                  color = "glom layer", opacity = 0.2,
                  text = ~text, hovertemplate = paste('%{text}'),
                  marker = list(size = 5, color = "grey")) %>%
        add_trace(data = plot_p50,
                  x = ~AntPos, y = ~MedLat, z = ~VenDor,
                  color = ~p50_rank,
                  text = ~text, hovertemplate = paste('%{text}'),
                  marker = list(size = 6, reversescale = T,
                                line = list(color = 'black', width = 0.5))) %>%
        layout(scene = list(xaxis = list(title = '(1)Anterior-Posterior(23)'),
                            yaxis = list(title = '(1)Medial-Lateral(22)'),
                            zaxis = list(title = '(22)Dorsal-Ventral(1)')),
               title = paste(rankx, "top ranked posterior median voxels for", gene))

    return(rank_plot)
} #end RankPosterior3D

ui <- fluidPage(
    titlePanel("OR-OB transcriptional mapping probabilities"),
    sidebarLayout(
        sidebarPanel(
            autocomplete_input("userOR",
                               "OR name",
                               options = gene_names,
                               max_options = 100),
            sliderInput("userRanks",
                        "Ranked probability voxels to show",
                        min = 1,
                        max = 2011,
                        value = 2011)
        ), mainPanel(
            plotlyOutput("OB3D")
        )
    )
) #end ui

server <- function(input, output) {
    output$OB3D <- renderPlotly(RankPosterior3D(gene = input$userOR,
                                              rankx = input$userRanks))
}

#Run the application
shinyApp(ui = ui, server = server)
