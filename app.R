# OBmap Shiny App
# Kevin Zhu
# Matsunami Lab, Duke University

# define UI for application that plots 3D voxels with OR glomeruli probabilities and assigned positions

#v2, 2022-02-10
# a. view additional OR information such as aliases, class, OE zone, etc.
# b. view predictions for OR of interest

#v1, 2021-09-15
# allows user to
# 1. select an OR of interest and number of top probability voxels to view interactively

#load packages
library(shiny)
library(dqshiny)
library(plotly)
library(dplyr)
library(readr)

#load model output and get gene names
probabilities <- read_csv("ob3d_posteriors_210808v1.csv") #71MB
assignments <- read_csv("ob3d_assignments_210808v1.csv") #4MB, includes all assignments, even those under p0005
info <- read_csv("info_210808v1.csv")

gene_names <- probabilities %>%
    select(olfrname) %>%
    arrange(olfrname) %>%
    unique() %>%
    pull(olfrname)

#given a gene and a cutoff number of top ranked voxels,
#rank and filter probabilities then plot interactively in 3D
RankGlom3D <- function(gene, rankx = 2011) {
    plot_p50 <- probabilities %>%
        filter(olfrname == gene) %>%
        mutate(Probability_Rank = min_rank(desc(p50))) %>%
        filter(Probability_Rank <= rankx) %>%
        mutate(text = paste('Rank:', Probability_Rank,
                            '<br>p50:', signif(p50, digits = 3),
                            '<br>AP:', AntPos,
                            '<br>VD:', VenDor,
                            '<br>ML:', MedLat))

    plot_glom <- assignments %>%
        filter(olfrname == gene) %>%
        mutate(Probability_Rank = min_rank(desc(p50)),
               text = paste('Rank:', Probability_Rank,
                            '<br>p50:', signif(p50, digits = 3),
                            '<br>AP:', AntPos,
                            '<br>VD:', VenDor,
                            '<br>ML:', MedLat),
               gsize = case_when(p50 == clustmaxp ~ 8,
                                 p50 != clustmaxp ~ 6))

    blank_p50 <- probabilities %>%
        filter(olfrname == gene) %>%
        mutate(Probability_Rank = min_rank(desc(p50))) %>%
        mutate(text = paste('Rank:', Probability_Rank,
                            '<br>p50:', signif(p50, digits = 3),
                            '<br>AP:', AntPos,
                            '<br>VD:', VenDor,
                            '<br>ML:', MedLat))

    rank_plot <- plot_ly(type = "scatter3d", mode = "markers") %>%
        add_trace(name = "OB shell",
                  data = blank_p50,
                  x = ~AntPos, y = ~MedLat, z = ~VenDor,
                  color = "OB shell", opacity = 0.2,
                  text = ~text, hovertemplate = paste('%{text}'),
                  marker = list(size = 5, color = "grey")) %>%
        add_trace(name = "Probabilities",
                  data = plot_p50,
                  x = ~AntPos, y = ~MedLat, z = ~VenDor,
                  color = ~Probability_Rank,
                  text = ~text, hovertemplate = paste('%{text}'),
                  marker = list(size = 6, reversescale = T,
                                line = list(color = 'black', width = 0.5))) %>%
        add_trace(name = "Assigned Position",
                  data = plot_glom,
                  x = ~AntPos, y = ~MedLat, z = ~VenDor,
                  color = "Assigned Position",
                  text = ~text, hovertemplate = paste('%{text}'),
                  marker = list(size = 8, color = "red",
                                line = list(color = 'black', width = 0.5))) %>%
        layout(scene = list(xaxis = list(title = '(1)Anterior-Posterior(23)'),
                            yaxis = list(title = '(1)Medial-Lateral(22)'),
                            zaxis = list(title = '(22)Dorsal-Ventral(1)')),
               title = paste("Top ranked posterior median voxels for", gene))

    return(rank_plot)
} #end RankGlom3D

DisplayInfo <- function(gene) {
    class <- paste(strong("Class:"), info %>%
                       filter(olfrname == gene) %>%
                       pull(class))
    tanzone <- paste(strong("Tan 2018 OE Zone:"), info %>%
                         filter(olfrname == gene) %>%
                         pull(tz_val))
    momzone <- paste(strong("Zapiec 2020 OE Zone:"), info %>%
                         filter(olfrname == gene) %>%
                         pull(Zolfr_Momb) %>%
                         round(digits = 2))
    fisur <- paste(strong("FI Surface Enriched:"), info %>%
                         filter(olfrname == gene) %>%
                         pull(fisurface))
    apmp <- paste(strong("AP Mean Position:"), info %>%
                         filter(olfrname == gene) %>%
                         pull(ap_wavg) %>%
                      round(digits = 2))
    dvmp <- paste(strong("DV Mean Position:"), info %>%
                         filter(olfrname == gene) %>%
                         pull(dv_wavg) %>%
                      round(digits = 2))
    mlmp <- paste(strong("ML Mean Position:"), info %>%
                         filter(olfrname == gene) %>%
                         pull(ml_wavg) %>%
                      round(digits = 2))
    paste(strong(gene), class, tanzone, momzone,
          fisur, apmp, dvmp, mlmp, sep = "<br>")
} #end DisplayInfo


ui <- fluidPage(
    titlePanel("OR-OB transcriptional mapping"),
    h5("Decoding the olfactory map: targeted transcriptomics link olfactory receptors to glomeruli"),
    h5("Zhu et al. 2022"),
    sidebarLayout(
        sidebarPanel(
            p("Enter Olfr style name for mouse OR genes:"),
            autocomplete_input("userOR",
                               "OR name",
                               options = gene_names,
                               max_options = 100),
            sliderInput("userRanks",
                        "Ranked probability voxels to show",
                        min = 1,
                        max = 2011,
                        value = 2011),
            hr(style = "border-top: 1px solid #000000;"),
            htmlOutput("info")

        ), mainPanel(
            plotlyOutput("OB3D")
        )
    )
) #end ui

server <- function(input, output) {
    output$OB3D <- renderPlotly(RankGlom3D(gene = input$userOR,
                                              rankx = input$userRanks))

    output$info <- renderText(DisplayInfo(gene = input$userOR))
} #end server

#Run the application
shinyApp(ui = ui, server = server)