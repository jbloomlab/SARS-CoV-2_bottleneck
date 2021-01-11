#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

## ==== Install Required Packages ==== ##

## List of all packages needed -- non-BioManager installed
library(shiny)
library(tidyverse)

## ==== Read in Pileup data ==== ##

combined.pileup.df = readRDS("combined-pileup-df.RDS")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Interactive Variant Calling"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("AF",
                        "Minimum Allele Frequency",
                        min = 0,
                        max = 0.05,
                        value = 0.001),
            
            numericInput(inputId = "DP",
                         "Depth", value = 10)
        ),
    
    

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("SNPs.Replicates")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$SNPs.Replicates <- renderPlot({
        
        ## ======= ##
        AF = input$AF
        DP = input$DP
        
        combined.pileup.df %>% 
            filter(Percent == "100") %>% 
            ggplot(aes(x = log10(AF.one), y = log10(AF.two))) + 
            facet_wrap(~Accession, ncol = 4) +
            geom_point(size = 3, alpha = 0.5, col = "#9c9c9c") +
            geom_point(data = filter(combined.pileup.df, DP.one >= DP & DP.two >= DP, AF.one >= AF & AF.two >= AF), 
                       mapping = aes(x = log10(AF.one), y = log10(AF.two)), size = 3, col = "#a30000") +
            geom_smooth(data = filter(combined.pileup.df, DP.one >= DP & DP.two >= DP, AF.one >= AF & AF.two >= AF), 
                        mapping = aes(x = log10(AF.one), y = log10(AF.two)),
                        method='lm', se = F, col = "#a30000") +
            geom_hline(yintercept = log10(AF), size = 0.5, linetype = 2) +
            geom_vline(xintercept = log10(AF), size = 0.5, linetype = 2) +
            xlab("Allele Frequency Replicate One") +
            ylab("Allele Frequency Replicate Two") +
            theme_classic() +
            theme(legend.position="bottom", legend.box = "horizontal") +
            theme(legend.box.background = element_rect(colour = "black")) +
            theme(text=element_text(size=18,  family="Helvetica")) +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(panel.grid.major.y = element_line(colour="grey", linetype="dashed")) +
            theme(panel.background = element_rect(fill = NA, color = "black"))
        
        ## ====== ##
        
    }, width = 900, height = 900)
}

# Run the application 
shinyApp(ui = ui, server = server)
