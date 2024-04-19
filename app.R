options(shiny.reactlog = TRUE)

library(rhandsontable)
library(shiny)
library(DT)
library(ggplot2)
library(magrittr)
library(data.table)
library(rhandsontable)
library(dplyr)
library(tidyverse)

life_table = data.frame(Age = c(1), Survival = c(0), Fecundity = c(0))
ui <- fluidPage(
  titlePanel("EEA Course Exercise"),
  tabsetPanel(id="out_disp",
              tabPanel("Life Table",
                       sidebarLayout(
                         sidebarPanel(
                           numericInput("rownumb", "Number of age classes", 1, min = 0, max = 50),
                           actionButton("add_btn", "Add")
                         ),
                         mainPanel(rHandsontableOutput("shiny_table", width = 400),
                                   actionButton("calc_btn", "Calculate!"))
                       )),
              # tabPanel("Plot", plotOutput("plot1"))
              tabPanel("Complete Life Table", value = "tab2",
                       htmlOutput("r_value"),
                       tableOutput("view_12"),
                       downloadButton("downloadData", "Download!")
              ),
              tabPanel("Plots",
                       sidebarLayout(
                         sidebarPanel(
                           uiOutput("var_choices")),
                         mainPanel(plotOutput("plot2", width = "75%")
                                   ,
                                   downloadButton("downloadplots", "Download!")
                         )
                       )
              ),
              tabPanel(
                "Hamilton", plotOutput("plot3", width = "75%")
                ,
                downloadButton("downloadhamil", "Download!")
              ),
              tabPanel(
                "Selection", plotOutput("plot1", width = "75%")
                ,
                downloadButton("downloadselgrad", "Download!")
              )
              # ,
              # tabPanel("Leslie Matrix", tableOutput("leslie_mat"))
  )
)


server <- function(input, output, session) {
  
  life_table <- reactiveVal(life_table)
  
  
  observeEvent(input$add_btn, {
    chkrow <- input$rownumb 
    # t1 <- data.frame(matrix(data=0, nrow=input$rownumb, ncol=3))
    t1 <- data.frame(matrix(data=0, nrow=chkrow, ncol=3))
    colnames(t1) <- c("Age", "Survival", "Fecundity")
    # t1$Age <- seq(1,input$rownumb, 1)
    t1$Age <- seq(1,chkrow, 1)
    t1[nrow(t1),2] <- NA
    t = tail(rbind(life_table(), t1),chkrow)
    life_table(t)
  })
  
  
  output$shiny_table <- renderRHandsontable({
    rhandsontable(life_table(), selection = 'multiple', digits=3, width=400, options = list(dom = 't', scrollY='400px', "pageLength" = 60), rowHeaders = NULL, editable = T)
  })
  
  observeEvent(input$calc_btn, {
    globals <-  hot_to_r(input$shiny_table) 
    NRows <- nrow(globals)
    r_range<- c(-1, 1)
    add_globals <- reactiveValues(df=globals)
    for(i in 1:NRows) {
      add_globals$df$qx[i] <- 1 - add_globals$df$Survival[i]
      add_globals$df$mux[i] <- -log(add_globals$df$Survival[i])}
    eul_lot <- reactiveVal(value = 0.1)
    p <- head(add_globals$df$Survival,NRows-1)
    m <- add_globals$df$Fecundity
    A <- diag(p)
    A <- cbind(A, rep(0,NRows-1))
    A <- rbind(m, A)
    r_val <- round(log(Re(eigen(A)$values[1])),3)
    eul_lot(r_val)
    
    new_tbl <- reactiveValues(dat=add_globals$df)
    sel_grad <- function(i){
      for (i in 1:NRows) {
        new_tbl$dat$lxmx[i] <- (new_tbl$dat$lx[i])*(new_tbl$dat$Fecundity[i])
        new_tbl$dat$erx[i] <- exp(-(eul_lot())*(new_tbl$dat$Age[i]))
        new_tbl$dat$product_erx_lxmx[i] <- (new_tbl$dat$lxmx[i])*(new_tbl$dat$erx[i])
        new_tbl$dat$Tgen[i] <- (new_tbl$dat$Age[i])*(new_tbl$dat$product_erx_lxmx[i])
        new_tbl$dat$T_sum <- sum((new_tbl$dat$Tgen), na.rm = TRUE)
        T_sum <- sum((new_tbl$dat$Tgen), na.rm = TRUE)
        for (i in 1:NRows){
          new_tbl$dat$s_mu_x[i] <- sum((new_tbl$dat$product_erx_lxmx[i+1:NRows]), na.rm=TRUE)/T_sum #from Hamilton's paper
          new_tbl$dat$s_mx[i] <- ((new_tbl$dat$lx[i])*(new_tbl$dat$erx[i]))/T_sum
          new_tbl$dat$s_mu_x_sel[i] <- sum((new_tbl$dat$product_erx_lxmx[i+1:NRows]), na.rm=TRUE) #from Hamilton's paper
          new_tbl$dat$s_mx_sel[i] <- (new_tbl$dat$lx[i])*(new_tbl$dat$erx[i])
        }
      }
      return(new_tbl$dat)
    }
    last_df <-  new_tbl$dat %>% sel_grad()
    last_df1 <- last_df %>% select(1:6, 11:15)
    colnames(last_df1) <- c("Age", "Survival","Fecundity", "Cumulative Survival (lx)", "Death Rate (qx)", "Mortality Rate (mu_x)", "Generation Time T", "Hamilton's Mortality Gradient", "Hamilton's Fecundity Gradient", "Mortality Selection", "Fecundity Selection")
    last_df1$r_value <- eul_lot()
    last_df1 <- last_df1 %>% mutate_if(is.numeric, round, digits=3)
    names(last_df)[12] <- "HamMortalitySel"
    names(last_df)[13] <- "HamFecunditySel"
    names(last_df)[14] <- "MortalitySel"
    names(last_df)[15] <- "FecunditySel"
    long <- reshape2::melt(last_df, id.vars = c("Age"))
    long1 <- subset(long, variable == "HamMortalitySel" | variable == "HamFecunditySel")
    long2 <- subset(long, variable == "MortalitySel" | variable == "FecunditySel")
    
    output$r_value <- renderText({paste("The value of r for which Euler-Lotka equation is satisfied is", "<b>", eul_lot(), "</b>")})
    output$view_12 <- renderTable({last_df1}, digits=3 )
    updateTabsetPanel(session, "out_disp",
                      selected = "tab2")
    fig_grad <- ggplot(long1, aes(Age, value, col=variable)) + geom_line(size=1.5,position=position_dodge(width=0.2)) + theme_bw() + theme(legend.text=element_text(size=16), axis.title=element_text(size=16), axis.text=element_text(size=12), legend.title = element_blank()) + labs(x = "Age", y = "Hamilton's Gradients") + scale_color_discrete(labels=c("Mortality", "Fecundity"))
    fig_sel <- ggplot(long2, aes(Age, value, col=variable)) + geom_line(size=1.5,position=position_dodge(width=0.2)) + theme_bw() + theme(legend.text=element_text(size=16), axis.title=element_text(size=16), axis.text=element_text(size=12), legend.title = element_blank()) + labs(x = "Age", y = "Selection Gradients") + scale_color_discrete(labels=c("Mortality", "Fecundity"))
    
    output$plot3 <- renderPlot({
      # plot1 <- ggplot(long1, aes(Age, value, col=variable)) + geom_line(size=1.5,position=position_dodge(width=0.2)) + theme_bw() + theme(legend.text=element_text(size=16), axis.title=element_text(size=16), axis.text=element_text(size=12)) + labs(x = "Age", y = "Selection") + scale_color_discrete(labels=c("Mortality", "Fecundity"))
      # plot1
      print(fig_grad)
    })
    
    output$plot1 <- renderPlot({
      # plot1 <- ggplot(long1, aes(Age, value, col=variable)) + geom_line(size=1.5,position=position_dodge(width=0.2)) + theme_bw() + theme(legend.text=element_text(size=16), axis.title=element_text(size=16), axis.text=element_text(size=12)) + labs(x = "Age", y = "Selection") + scale_color_discrete(labels=c("Mortality", "Fecundity"))
      # plot1
      print(fig_sel)
    })
    
    choices <- list("Age", "Survival", "Fecundity", "Cumulative Survival (lx)", "Death Rate (qx)", "Mortality Rate (mu_x)")
    output$var_choices <- renderUI({
      tagList(
        selectInput('xcol', 'X Variable', choices),
        selectInput('ycol', 'Y Variable', choices)
      )
    })
    fin_df <- reactive({last_df1[, c(input$xcol, input$ycol)]})
    fig_plots <-  reactive({plot(fin_df(), type = "l")})
    output$plot2 <- renderPlot({
      # par(mar = c(4, 4, 4, 4))
      # plot(fin_df(), type = "l")
      print(plot(fin_df(), type = "l"))
    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("EEA_TotalLifeTable", ".csv", sep = "")
      },
      content = function(file) {
        write.csv(last_df1, file, row.names = FALSE)
      }
    )
    output$downloadhamil <- downloadHandler(
      filename = function() {
        paste("EEA_HamiltonPlots", ".png", sep = "")
      },
      content = function(file) {
        # png(file)
        # print(fig_grad())
        # dev.off()
        ggsave(file, plot = fig_grad, device = "png", width = 7, height = 5, dpi = 300, units = "in")
      }
    )
    output$downloadselgrad <- downloadHandler(
      filename = function() {
        paste("EEA_SelectionPlots", ".png", sep = "")
      },
      content = function(file) {
        # png(file)
        # print(fig_grad())
        # dev.off()
        ggsave(file, plot = fig_sel, device = "png", width = 7, height = 5, dpi = 300, units = "in")
      }
    )
    output$downloadplots <- downloadHandler(
      filename = function() {
        paste("EEA_", input$xcol, "Vs", input$ycol, "Plot.png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = fig_plots(), device = "png", width = 6, height = 6, dpi = 300, units = "in")
        # png(file)
        # print(fig_plots())
        # dev.off()  
      }
    )
  })
}

shinyApp(ui = ui, server = server)   
