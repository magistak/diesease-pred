library(shiny)
library(tidyverse)
library(janitor)
library(forcats)
library(purrr)

met_min_tbl <- read_rds("data/met_min_tbl_v2.rds")
v_hr_tbl <- read_rds("data/v_hr_tbl.rds")
nmr_vector <- c("Albumin","Creatinine","HDL Cholesterol","Clinical LDL Cholesterol","Apolipoprotein A1","Apolipoprotein B", "Total Cholesterol", "Total Triglycerides", "non_HDL") |> make_clean_names()
clin_vector <- c("Albumin_clinical", "Creatinine_clinical", "HDL cholesterol_clinical", "LDL direct_clinical", "Apolipoprotein A_clinical", "Apolipoprotein B_clinical", "Cholesterol_clinical", "Triglycerides_clinical", "non_HDL_clinical") |> make_clean_names()
disease_vector <- c("all-cause mortality", "T2 Diabetes", "COPD", "Liver Disease", "cvd") |> make_clean_names()
biom_vector <- c(nmr_vector, clin_vector)
filters_tbl <- read_rds("data/filters_tbl.rds")
# add percents significant


heat_min <- 

heat_tbl <- 
met_min_tbl |> 
  left_join(v_hr_tbl) |> 
  mutate(p_signif = lrt_p<0.05,
         v_shape = hr_left < 1 & hr_right > 1,
         only_signif=ifelse(p_signif,lrt_p,NA),
         vp = ifelse(v_shape, lrt_p, NA))

ui <- fluidPage(
  titlePanel(
    h1(textOutput("title_panel"), align = "center")
  ),
  sidebarLayout(
    sidebarPanel(width = 2,
                 # Add your sidebar widgets here
                 selectInput("sex_filter_i", "Sex Filter:", choices = c("both", "men", "women")),
                 selectInput("health_filter_i", "Health Filter:", choices = c("is_healthy_all", "is_healthy_7", "is_healthy_3", "is_healthy_1")),
                 textOutput("no_signif"),
                 textOutput("cohort_size")
    ),
    mainPanel(
      fluidRow(
        column(
          width = 6,
          # Add the first clickable heatmap
          plotOutput("heatmap_nmr", click = "heatmap1"),
          # Add the second clickable heatmap
          plotOutput("heatmap_clin", click = "heatmap2")
        ),
        column(
          width = 3,
          # Display the clicked cell values
          imageOutput("pic_spline_nmr"),
          imageOutput("pic_spline_clin")
        ),
        column(
          width = 3,
          imageOutput("pic_v_nmr"),
          imageOutput("pic_v_clin")
        )
      )
    )
  )
)

server <- function(input, output){
  
  heatmap_data <- reactive({

    heat_tbl |> 
      filter(sex_filter == input$sex_filter_i,
             health_filter == input$health_filter_i) |> 
      arrange(match(disease, disease_vector), match(metabolite, biom_vector))
  })
  
  heatmap_plot1 <- reactive({
    
    heatmap_data() |>
      filter(nmr_clin == "nmr") |>
      ggplot(aes(x = fct_inorder(metabolite), y = fct_inorder(disease))) +
      geom_tile(aes(fill = vp, color = lrt_p < 0.05 & !is.na(lrt_p)), linewidth = 1.5) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        title = sprintf("Sex: %s\nHealth filter: %s\nHeatmap: nmr", 
                        input$sex_filter_i, input$health_filter_i),
        x = "",
        y = ""
      ) +
   #   scale_fill_manual(values = c("above" = "green", "below" = "red", "other" = "blueviolet")) +
      scale_fill_gradient(low = "blue", high = "white", na.value = "gray") +
      theme(text = element_text(size = 20, face="bold")) +
      guides(fill = guide_colorbar(
        ticks = FALSE,
        frame.colour = "black",
        title = "p"
      )) +
      guides(color=guide_legend(title="p < 0.05"))
  })
  
  heatmap_plot2 <- reactive({
    
    heatmap_data() |>
      filter(nmr_clin == "clinical") |>
      ggplot(aes(x = fct_inorder(metabolite), y = fct_inorder(disease))) +
      geom_tile(aes(fill = vp, color = lrt_p < 0.05 & !is.na(lrt_p)), linewidth = 1.5) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      labs(
        title = sprintf("Sex: %s\nHealth filter: %s\nHeatmap: clinical", 
                        input$sex_filter_i, input$health_filter_i),
        x = "",
        y = ""
      ) +
      #   scale_fill_manual(values = c("above" = "green", "below" = "red", "other" = "blueviolet")) +
      scale_fill_gradient(low = "blue", high = "white", na.value = "gray") +
      theme(text = element_text(size = 20, face="bold")) +
      guides(fill = guide_colorbar(
        ticks = FALSE,
        frame.colour = "black",
        title = "p"
      )) +
      guides(color=guide_legend(title="p < 0.05"))
  })
  
  # Output the heatmap_nmr plot
  output$heatmap_nmr <- renderPlot({
    heatmap_plot1()
  })
  
  output$heatmap_clin <- renderPlot({
    heatmap_plot2()
  })
  
  clicked_cell <- reactiveValues(
    disease = 1,
    metabolite = 1
  )
  
  pics <- reactiveValues(
    sex_filter = "both",
    health_filter = "is_healthy_all",
    disease_index = 1,
    met_index = 1
  )
  
  observeEvent(input$heatmap1, {
    clicked_cell$disease <- floor(input$heatmap1$y-0.5)+1
    clicked_cell$metabolite <- floor(input$heatmap1$x-0.5)+1
    pics$disease_index <- floor(input$heatmap1$y-0.5)+1
    pics$met_index <- floor(input$heatmap1$x-0.5)+1
  })
  
  observeEvent(input$heatmap2, {
    clicked_cell$disease <- floor(input$heatmap2$y-0.5)+1
    clicked_cell$metabolite <- floor(input$heatmap2$x-0.5)+1
    pics$disease_index <- floor(input$heatmap2$y-0.5)+1
    pics$met_index <- floor(input$heatmap2$x-0.5)+1
  })
  
  output$title_panel <- renderText({
    paste(nmr_vector[clicked_cell$metabolite], disease_vector[clicked_cell$disease])
  })
  
  
  observeEvent(input$sex_filter_i, {
    pics$sex_filter <- input$sex_filter_i
  })
  
  observeEvent(input$health_filter_i, {
    pics$health_filter <- input$health_filter_i
  })
  
  output$pic_spline_nmr <- renderImage({
    file_path <- sprintf("%s%s#%s#%s#%s#%s#%s.png","data/pics/spline/",
                         disease_vector[pics$disease_index], nmr_vector[pics$met_index],
                         pics$sex_filter, pics$health_filter, "no", "nmr")
    list(src = file_path, contentType = "image/png", alt = "PNG Plot",width = "120%", height = "100%")
  }, deleteFile = FALSE)
  
 
  output$pic_spline_clin <- renderImage({
    file_path <- sprintf("%s%s#%s#%s#%s#%s#%s.png","data/pics/spline/",
                         disease_vector[pics$disease_index], clin_vector[pics$met_index],
                         pics$sex_filter, pics$health_filter, "no", "clin")
    list(src = file_path, contentType = "image/png", alt = "PNG Plot",width = "120%", height = "100%")
  }, deleteFile = FALSE)
  
  output$pic_v_nmr <- renderImage({
    file_path <- sprintf("%s%s#%s#%s#%s#%s#%s.png","data/pics/v_shape/",
                         disease_vector[pics$disease_index], nmr_vector[pics$met_index],
                         pics$sex_filter, pics$health_filter, "no", "nmr")
    list(src = file_path, contentType = "image/png", alt = "PNG Plot",width = "120%", height = "100%")
  }, deleteFile = FALSE)
  
  
  output$pic_v_clin <- renderImage({
    file_path <- sprintf("%s%s#%s#%s#%s#%s#%s.png","data/pics/v_shape/",
                         disease_vector[pics$disease_index], clin_vector[pics$met_index],
                         pics$sex_filter, pics$health_filter, "no", "clinical")
    list(src = file_path, contentType = "image/png", alt = "PNG Plot",width = "120%", height = "100%")
  }, deleteFile = FALSE)
  
  output$no_signif <- renderText({
  nrh <- nrow(heatmap_data() |> 
    filter(nmr_clin == "nmr"))
  signif_nrh <- nrow(heatmap_data() |> 
                       filter(!is.na(only_signif)) |> 
                       filter(nmr_clin == "nmr"))

    sprintf("Percent of significants nmr: %s",signif(signif_nrh/nrh*100,4))
  })
  
  output$cohort_size <- renderText({
    a_s <- if (pics$sex_filter == "both") c("0", "1") else if (pics$sex_filter == "men") 1 else 0
    
    eid_filter <- 
      filters_tbl |> 
      filter(!!sym(pics$health_filter) == 1) |> 
      filter(sex %in% a_s)
    
    sprintf("Cohort size: %s",length(eid_filter$eid))

      })
  
}

shinyApp(ui, server)
