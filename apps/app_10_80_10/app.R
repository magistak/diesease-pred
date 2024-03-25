library(shiny)
library(tidyverse)
library(forcats)
library(janitor)

heat_tbl <- read_rds("data/heat_tbl.rds")
nmr_vector <- c("Albumin","Creatinine","HDL Cholesterol","Clinical LDL Cholesterol","Apolipoprotein A1","Apolipoprotein B", "Total Cholesterol", "Total Triglycerides", "non_HDL") |> make_clean_names()
clin_vector <- c("Albumin_clinical", "Creatinine_clinical", "HDL cholesterol_clinical", "LDL direct_clinical", "Apolipoprotein A_clinical", "Apolipoprotein B_clinical", "Cholesterol_clinical", "Triglycerides_clinical", "non_HDL_clinical") |> make_clean_names()
disease_vector <- c("all-cause mortality", "T2 Diabetes", "COPD", "Liver Disease", "cvd") |> make_clean_names()


ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(width = 2,
      # Add your sidebar widgets here
      selectInput("sex_filter_i", "Sex Filter:", choices = c("men", "women", "both")),
      selectInput("health_filter_i", "Health Filter:", choices = c("is_healthy_1", "is_healthy_3", "is_healthy_7", "is_healthy_all")),
      selectInput("bmi_filter_i", "BMI Filter:", choices = c("yes", "no"))
    ),
    mainPanel(
      fluidRow(
        column(
          width = 8,
          # Add the first clickable heatmap
          plotOutput("heatmap_clin", click = "heatmap1"),
          # Add the second clickable heatmap
          plotOutput("heatmap_nmr", click = "heatmap2")
        ),
        column(
          width = 4,
          # Display the clicked cell values
          imageOutput("plot", height = "10%", width = "10%")
        )
      )
    )
  )
)

server <- function(input, output) {
  # Reactive expression for the heatmap data
  heatmap_data <- reactive({
    sex_filter_var <- input$sex_filter_i
    health_filter_var <- input$health_filter_i
    bmi_filter_var <- input$bmi_filter_i

      heat_tbl |> 
      filter(sex_filter == sex_filter_var,
             health_filter == health_filter_var,
             bmi_filter == bmi_filter_var) |> 
      arrange(match(disease, disease_list), match(metabolite, biom_list))
  })
  
  # Reactive expression for the heatmap_clin plot
  heatmap_clin_plot <- reactive({
    
    heatmap_data() |>
      filter(nmr_clin == "clinical") |>
      ggplot(aes(x = fct_inorder(metabolite), y = fct_inorder(disease), fill = heat_group)) +
      geom_tile(color = "black") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        title = sprintf("Sex: %s\nHealth filter: %s\nBmi filter: %s\nHeatmap: Clin", 
                        sex_filter_var, health_filter_var, bmi_filter_var),
        x = "",
        y = ""
      ) +
      scale_fill_manual(values = c("above" = "green", "below" = "red", "other" = "blueviolet")) +
      theme(text = element_text(size = 10))
  })
  
  # Reactive expression for the heatmap_nmr plot
  heatmap_nmr_plot <- reactive({
    
    heatmap_data() |>
      filter(nmr_clin == "nmr") |>
      ggplot(aes(x = fct_inorder(metabolite), y = fct_inorder(disease), fill = heat_group)) +
      geom_tile(color = "black") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        title = sprintf("Sex: %s\nHealth filter: %s\nBmi filter: %s\nHeatmap: NMR", 
                        input$sex_filter_i, input$health_filter_i, input$bmi_filter_i),
        x = "",
        y = ""
      ) +
      scale_fill_manual(values = c("above" = "green", "below" = "red", "other" = "blueviolet")) +
      theme(text = element_text(size = 10))
  })
  
  # Output the heatmap_clin plot
  output$heatmap_clin <- renderPlot({
    heatmap_clin_plot()
  })
  
  # Output the heatmap_nmr plot
  output$heatmap_nmr <- renderPlot({
    heatmap_nmr_plot()
  })
  
  # Reactive values to store the clicked cell values
  clicked_cell <- reactiveValues(
    disease = NULL,
    metabolite = NULL,
    nmr_clin = ""
  )
  
  # Observe the click event on heatmap_clin
  observeEvent(input$clin_click, {
    clicked_cell$disease <- floor(as.numeric(input$clin_click$y)) 
    clicked_cell$metabolite <- floor(as.numeric(input$clin_click$x)) + 1
    clicked_cell$nmr_clin <- "clin"
  })
  
  # Observe the click event on heatmap_nmr
  observeEvent(input$nmr_click, {
    clicked_cell$disease <- floor(as.numeric(input$nmr_click$y)) 
    clicked_cell$metabolite <- floor(as.numeric(input$nmr_click$x)) + 1
    clicked_cell$nmr_clin <- "nmr"
  })
  
  # Create the PNG file path based on the clicked cell values
  png_file_path <- reactive({
    req(clicked_cell$disease, clicked_cell$metabolite, clicked_cell$nmr_clin)
    
    if (clicked_cell$nmr_clin == "clin") {
      metabolite_name <- clin_vector[clicked_cell$metabolite]
    } else {
      metabolite_name <- nmr_vector[clicked_cell$metabolite]
    }
    
    disease_name <- disease_vector[clicked_cell$disease]
    
    sex_filter_var <- input$sex_filter_i
    health_filter_var <- input$health_filter_i
    bmi_filter_var <- if (input$bmi_filter_i == "yes") "yes" else "no"
    
    file_path <- sprintf("%s%s#%s#%s#%s#%s#%s.png","data/spline_pics/",
                         disease_name, metabolite_name,
                         sex_filter_var, health_filter_var, bmi_filter_var,
                         clicked_cell$nmr_clin)
    
    return(file_path)
  })
  
  #  Output the PNG plot
  output$plot <- renderImage({
    req(clicked_cell$disease, clicked_cell$metabolite, clicked_cell$nmr_clin)
    list(src = png_file_path(), contentType = "image/png", alt = "PNG Plot")
  }, deleteFile = FALSE)
  
  # Output the heatmap_clin plot
  output$heatmap_clin <- renderPlot({
    heatmap_clin_plot()
  })
  
  # Output the heatmap_nmr plot
  output$heatmap_nmr <- renderPlot({
    heatmap_nmr_plot()
  })


}

shinyApp(ui = ui, server = server)