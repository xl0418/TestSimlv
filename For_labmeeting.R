library(shiny)
library(plotly)



death_rate <- function(m, temp, temp_opt, sigma){
  return(m + (temp_opt - temp)^2 / sigma^2)
}


sim_rc_varyingm <- function(tt, mu1, mu2,
                   v1, v2,
                   r1in,
                   N10, N20, R10,
                   m,
                   temp, temp_opt1, temp_opt2, sigma
) {
  N1 <- c(N10)
  N2 <- c(N20)
  R1 <- c(R10)
  m1 <- death_rate(m, temp, temp_opt1, sigma)
  m2 <- death_rate(m, temp, temp_opt2, sigma)
  for(ti in 1:tt){
    temp_N1 <- N1[ti]
    temp_N2 <- N2[ti]
    temp_R1 <- R1[ti]
    for(tii in 1:100){
      dN1 <- (mu1 * temp_R1 - m1) * temp_N1 / 100
      dN2 <- (mu2 * temp_R1 - m2) * temp_N2 / 100
      dR1 <- ((r1in - temp_R1) - v1 * temp_N1 - v2 * temp_N2) / 100

      temp_N1 <- max(temp_N1 + dN1, 0)
      temp_N2 <- max(temp_N2 + dN2, 0)
      temp_R1 <- max(temp_R1 + dR1, 0)

    }

    N1 <- c(N1, temp_N1)
    N2 <- c(N2, temp_N2)

    R1 <- c(R1, temp_R1)
  }
  return(data.frame(N1 = N1, N2 = N2, R1 = R1, cc = 1:(tt+1)))

}


growth_rate <- function(mu, temp, temp_opt, sigma){
  return(mu * exp( -(temp_opt - temp)^2 / sigma^2))
}


sim_rc_varyingb <- function(tt, mu1, mu2,
                   v1, v2,
                   r1in,
                   N10, N20, R10,
                   m,
                   temp, temp_opt1, temp_opt2, sigma, tau
) {
  N1 <- c(N10)
  N2 <- c(N20)
  R1 <- c(R10)

  for(ti in 1:tt){
    temp_N1 <- N1[ti]
    temp_N2 <- N2[ti]
    temp_R1 <- R1[ti]
    for(tii in 1:100){
      temp <- sin(2 * pi * (ti + tii / 100) / tau)
      mu1_g <- growth_rate(mu1, temp, temp_opt1, sigma)
      mu2_g <- growth_rate(mu2, temp, temp_opt2, sigma)
      dN1 <- (mu1_g * temp_R1 - m) * temp_N1 / 100
      dN2 <- (mu2_g * temp_R1 - m) * temp_N2 / 100
      dR1 <- ((r1in - temp_R1) - v1 * mu1_g * temp_N1 * temp_R1 - v2 *mu2_g  * temp_N2 * temp_R1) / 100

      temp_N1 <- max(temp_N1 + dN1, 0)
      temp_N2 <- max(temp_N2 + dN2, 0)
      temp_R1 <- max(temp_R1 + dR1, 0)

    }

    N1 <- c(N1, temp_N1)
    N2 <- c(N2, temp_N2)

    R1 <- c(R1, temp_R1)
  }
  return(data.frame(N1 = N1, N2 = N2, R1 = R1, cc = 1:(tt+1)))

}


Sim_rc_ess <- function(tt, mu11, mu12, mu21, mu22,
                       q11, q12, q21, q22,
                       r1in, r2in,
                       N10, N20, R10, R20,
                       m) {
  N1 <- c(N10)
  N2 <- c(N20)
  R1 <- c(R10)
  R2 <- c(R20)
  for(ti in 1:tt){
    temp_N1 <- N1[ti]
    temp_N2 <- N2[ti]
    temp_R1 <- R1[ti]
    temp_R2 <- R2[ti]
    for(tii in 1:100){

      mu1 <- min(mu11 * temp_R1, mu12 * temp_R2)
      mu2 <- min(mu21 * temp_R1, mu22 * temp_R2)


      dN1 <- (mu1  - m) * temp_N1 / 100
      dN2 <- (mu2  - m) * temp_N2 / 100
      dR1 <- ((r1in - temp_R1) - q11 * mu1 * temp_N1 - q21 * mu2 * temp_N2) / 100
      dR2 <- ((r2in - temp_R2) - q12 * mu1 * temp_N1 - q22 * mu2 * temp_N2) / 100

      temp_N1 <- max(temp_N1 + dN1, 0)
      temp_N2 <- max(temp_N2 + dN2, 0)
      temp_R1 <- max(temp_R1 + dR1, 0)
      temp_R2 <- max(temp_R2 + dR2, 0)

    }

    N1 <- c(N1, temp_N1)
    N2 <- c(N2, temp_N2)

    R1 <- c(R1, temp_R1)
    R2 <- c(R2, temp_R2)
  }
  return(data.frame(N1 = N1, N2 = N2, R1 = R1, R2 = R2, cc = 1:(tt+1)))

}



# Define UI for app that draws a histogram ----
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  fluidRow(

    tabsetPanel(id = "tabset",
                tabPanel("Varying m",
                         fluidRow(
                           column(4,
                                  div(style="height: 80px;",sliderInput('tempopt1', HTML("T<sub>opt,1</sub>"), -0.2, min = -2., max = 2., step = 0.1)),
                                  div(style="height: 80px;",sliderInput('tempopt2', HTML("T<sub>opt,2</sub>"), 0.4, min = -2., max = 2., step = 0.1))
                           ),
                           column(4,
                                  div(style="height: 80px;", sliderInput('T', HTML("Temperature"), 0.5, min = -1.5, max = 1.5))
                           )
                         )
                ),
                tabPanel("Varying growth rates",
                         fluidRow(
                           column(4,
                                  div(style="height: 80px;",sliderInput('tempopt1', HTML("T<sub>opt,1</sub>"), -0.2, min = -2., max = 2., step = 0.1)),
                                  div(style="height: 80px;",sliderInput('tempopt2', HTML("T<sub>opt,2</sub>"), 0.4, min = -2., max = 2., step = 0.1))
                           ),
                           column(4,
                                  div(style="height: 80px;", sliderInput('tau', HTML("Period"), 3.5, min = 0.1, max = 20))
                           )
                         )
                )
    )),

  fluidRow(
    column(12,
           plotlyOutput("plot1", width=1000, height = 400),
           plotlyOutput("plot2", width=1000, height=400)
    ))

)

server <- function(input, output, session) {




    output$plot1 <- renderPlotly({
      if (input$tabset == "Varying m") {
        # sim lv model
        sim_result <- reactive({
          sim_rc_varyingm(100, 1.0, 1.0, 1.3, 1.3,
                          10, 0.1, 0.1, 1,
                          0.1, input$T, input$tempopt1, input$tempopt2, 1
          )
        })
        # browser()

        T <- seq(-1, 1, 0.01)

        m1 <-reactive({
          death_rate(m = 0.1, temp = T, temp_opt = input$tempopt1, sigma = 1)
        })

        m2 <-reactive({
          death_rate(m = 0.1, temp = T, temp_opt = input$tempopt2, sigma = 1)
        })


        gr <- data.frame(T = T, m1 = m1(), m2 = m2())
        # browser()
        p <- plot_ly(data = gr) %>%
          add_trace(x = ~T, y = ~m1, type = 'scatter', mode = 'lines', name = 'm1', line = list(color = 'rgb(205, 12, 24)', width = 4)) %>%
          add_trace(x = ~T, y = ~m2, type = 'scatter', mode = 'lines', name = 'm2', line = list(color = 'rgb(22, 96, 167)', width = 4))

        # add vertical line with input$T
        p <- p %>% add_trace(x = c(input$T, input$T), y = c(0, 2), type = 'scatter', mode = 'lines', name = 'Temperature',
                             line = list(color = 'black', width = 4, dash = 'dash'))

        # xlim
        # p <- p %>% layout(xaxis = list(range = c(0, 2)),
        #                   yaxis = list(range = c(0, 2)))
        p <- p %>% layout(xaxis = list(title = "R"),yaxis = list(title = "Death rate"))



        p
        }else if (input$tabset == "Varying growth rates") {
          # sim lv model
          sim_result <- reactive({
            sim_rc_varyingb(100, 1.0, 1.0, 1.3, 1.3,
                            10, 0.1, 0.1, 1,
                            0.1, 0, input$tempopt1, input$tempopt2, 1, input$tau
            )
          })
        }
        })


    output$plot2 <- renderPlotly({
      if (input$tabset == "Varying m") {
        # sim lv model
        sim_result <- reactive({
          sim_rc_varyingm(100, 1.0, 1.0, 1.3, 1.3,
                          10, 0.1, 0.1, 1,
                          0.1, input$T, input$tempopt1, input$tempopt2, 1
          )
        })
        # browser()
        p <- plot_ly(data = sim_result()) %>%
          add_trace(x = ~cc, y = ~N1, type = 'scatter', mode = 'lines', name = 'Species 1', line = list(color = 'rgb(205, 12, 24)', width = 4)) %>%
          add_trace(x = ~cc, y = ~N2, type = 'scatter', mode = 'lines', name = 'Species 2', line = list(color = 'rgb(22, 96, 167)', width = 4))


        # xlim
        # p <- p %>% layout(xaxis = list(range = c(0, 2)),
        #                   yaxis = list(range = c(0, 2)))

        p <- p %>% layout(xaxis = list(title = "Time"),yaxis = list(title = "Density"))


        p
      }

    })



}
# Create Shiny app ----
shinyApp(ui = ui, server = server)
