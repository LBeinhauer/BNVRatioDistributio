#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)


# Probability Density Functions of:

# uncorrelated central normal ratio #

# pdf of ratio Z = X/Y
p_uncor_cent <- function(z){
    (1/(pi))/(1 + z^2)
}




# uncorrelated non-central normal ratio

# pdf of ratio Z = X/Y, with non-zero means
p_uncor <- function(z, sig_y, sig_x, mu_x, mu_y){
    
    a <- sqrt((1/(sig_x^2))*(z^2) + (1/(sig_y^2)))
    
    b <- ((mu_x/(sig_x^2))*z + (mu_y/(sig_y^2)))
    
    c <- ((mu_x^2)/(sig_x^2)) + ((mu_y^2)/(sig_y^2))
    
    d <- exp(((b^2) - (c*(a^2)))/(2*(a^2)))
    
    p <- ((b*d)/(a^3)) * (1/(sqrt(2*pi)*sig_x*sig_y)) * (pnorm((b/a)) - pnorm(-(b/a))) + 
        (1/((a^2)*pi*sig_x*sig_y)) * exp(-(c/2))
    
    p
}


dbnvr <- function(z, sig_y = 1, sig_x = 1, mu_x = 0, mu_y = 0, cor = 0){
    
    # uncorrelated central
    if(cor == 0 & sum((c(mu_x, mu_y) == 0)) == 2){
        p <- p_uncor_cent(z)
        info <- "uncorrelated central"
    }
    
    # uncorrelated non-central
    if(cor == 0 & sum((c(mu_x, mu_y) == 0)) != 2){
        p <- p_uncor(z, sig_y, sig_x, mu_x, mu_y)
        info <- "uncorrelated non-central"
    }
    
    # correlated central
    if(cor != 0 & sum((c(mu_x, mu_y) == 0)) == 2){
        p <- p_cor_cent(z, sig_y, sig_x, cor)
        info <- "correlated central"
    }
    
    # correlated non-central
    if(cor != 0 & sum((c(mu_x, mu_y) == 0)) != 2){
        stop("PDF for Ratio of correlated, non-central bivariate normal distribution currently not implemented!")
    }
    
    c(p, info)
    
}





pbnvr <- function(z, sig_y = 1, sig_x = 1, mu_x = 0, mu_y = 0, cor = 0){
    
    # uncorrelated central
    if(cor == 0 & sum((c(mu_x, mu_y) == 0)) == 2){
        dens <- function(x){p_uncor_cent(x)}
        
        p <- integrate(dens, lower = -Inf, upper = z)
    }
    
    # uncorrelated non-central
    if(cor == 0 & sum((c(mu_x, mu_y) == 0)) != 2){
        dens <- function(x){p_uncor(x, sig_y, sig_x, mu_x, mu_y)}
        
        p <- integrate(dens, lower = -Inf, upper = z)
    }
    
    # correlated central
    if(cor != 0 & sum((c(mu_x, mu_y) == 0)) == 2){
        dens <- function(x){p_cor_cent(x, sig_y, sig_x, cor)}
        
        p <- integrate(dens, lower = -Inf, upper = z)
    }
    
    # correlated non-central
    if(cor != 0 & sum((c(mu_x, mu_y) == 0)) != 2){
        stop("PDF for Ratio of correlated, non-central bivariate normal distribution currently not implemented!")
    }
    
    p
    
}





# correlated central normal ratio

# pdf of ratio Z = X/Y, with zero means and correlated variables
p_cor_cent <- function(z, sig_y, sig_x, cor){
    
    alph <- cor*(sig_x/sig_y)
    
    bet <-  (sig_x/sig_y) * sqrt(1 - cor)
    
    p <- (1/pi) * (bet/(((z-alph)^2) + (bet^2)))
    
    p
    
}


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Density of Ratio of Bivariate Normal Variables X & Y"),

    # Sidebar with a slider input for number of bins 
    sidebarPanel(
            numericInput("mu_x",
                        "Population Mean(X):",
                        min = -100,
                        max = 100,
                        value = 10)
        ),
        
        sidebarPanel(
            numericInput("mu_y",
                        "Population Mean(Y):",
                        min = -100,
                        max = 100,
                        value = 10)
        ),
        
        sidebarPanel(
            numericInput("sig_x",
                         "Population SD(X):",
                         min = 0,
                         max = 30,
                         value = 3)
        ),
        
        sidebarPanel(
            numericInput("sig_y",
                         "Population SD(y):",
                         min = 0,
                         max = 30,
                         value = 3)
        ),
        
        sidebarPanel(
            numericInput("cor",
                         "Desired correlation of X & Y:",
                         min = -1,
                         max = 1,
                         value = 0)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )


# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        mat <- matrix(NA, nrow = length(seq(-3, 3, .1)), ncol = 2)
        
        for(i in 1:length(seq(-3, 3, .1))){
            z <- seq(-3, 3, .1)[i]
            mat[i,1] <- dbnvr(z = z, mu_x = input$mu_x, mu_y = input$mu_y, sig_x = input$sig_x, 
                              sig_y = input$sig_y, cor = input$cor)[1]
            mat[i,2] <- z
        }
        
        plot(mat[,2], mat[,1])
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
