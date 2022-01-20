setwd("C:/Users/Kenji/Desktop/rocker_proj/roc_out")

library(ggplot2)
library(data.table)
library(plotly)
library(reticulate)
library(shiny)
library(shinyalert)



#Locate the files I want and name them nicely
initial_data_gathering = function(directory){
  
#Find the aligned reads
df = data.table(path = list.files(path = directory, recursive = T, full.names = T, pattern = "aligned_reads.txt", ignore.case = F))
#postive or negative
df[, category := ifelse(grepl("/negative/", path), "Negative", "Positive")]
protein_names = unlist(strsplit(df$path, split = "/aligned_reads/"))[c(F,T)]
protein_names = unlist(strsplit(protein_names, split = "_aligned_reads.txt"))
df[, protein_name := protein_names]
rm(protein_names)

return(df)
}

#Use the file paths to load the files
load_reads_and_describe = function(rdf){
  
  #Load the aligned read files, parse annotations, assign file origin and concatenate
  read_data = rbindlist(lapply(rdf$path, function(x){
    
    protein_name = unlist(strsplit(x, split = "/aligned_reads/"))[c(F,T)]
    protein_name = unlist(strsplit(protein_name, split = "_aligned_reads.txt"))
    
    one_file = fread(x, header = F)
    #no hit files will return as nulls; this is not uncommon with poor negative targets
    
    if(nrow(one_file) > 0){
      
      #The only ones I care about are query, target, sstart, ssend, and bitscore
      colnames(one_file) = c("query", "tgt", "pid", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
      one_file[, pid := NULL]
      one_file[, length := NULL]
      one_file[, mismatch := NULL]
      one_file[, gapopen := NULL]
      one_file[, qstart := NULL]
      one_file[, qend := NULL]
      one_file[, evalue := NULL]
      
      annot_string = unlist(strsplit(one_file$query, split = ";"))
      
      one_file[, read_id        := annot_string[c(T, F, F, F, F, F)]]
      one_file[, origin_start   := annot_string[c(F, T, F, F, F, F)]]
      one_file[, origin_end     := annot_string[c(F, F, T, F, F, F)]]
      one_file[, strand         := annot_string[c(F, F, F, T, F, F)]]
      one_file[, origin_genome  := annot_string[c(F, F, F, F, T, F)]]
      one_file[, classifier     := annot_string[c(F, F, F, F, F, T)]]
     
      one_file[, origin := protein_name]
       
    }
    
    
  }))
  
  return(read_data)
  
}


#shiny code

#helper functions
{
  #OS-agnostic UI dir chooser
  choose_directory = function(caption = 'Select data directory') {
    if (exists('choose.dir')) {
      choose.dir(caption = caption)
    } else {
      easycsv::choose_dir()
    }
  }
  
  
  warning_plot = function(){
    warning_plot <- ggplot(data = NULL, aes(x = 1, y = 1, label = "There is no plot ready to load yet. You'll have to prepare a directory."))+
      geom_text(size = 6) +
      theme(panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
    return(warning_plot)
  }
  
  
}


rocker_ui <-function(){
  
  ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }
           "
      )
    )
  ),
  
  useShinyalert(),
  
  sidebarPanel(width = 3,
         #Selectable options go here
         
         actionButton('rocker_dir', 'Select a ROCker project', icon = icon("folder-open")),
         verbatimTextOutput("roc_dir_name"),
         
         actionButton('load_dir', "Load the ROCker directory", icon = icon("folder-upload")),
         
         checkboxGroupInput("proteins", label="Proteins to Show")
         
         ),
  mainPanel(id = "rocker_main",
            #Spacing
            fluidRow(
              column(12,
                     div(style = "height:60px;background-color: white;", "")
              )
            ),
            fluidRow(
              
              plotlyOutput("rocker_plot", height = "850px")
              
            )
            
            
  )
  
)
  
  return(ui)
  
}

# Server
rocker_server <- function(input, output, session) {
  directory = getwd()
  output$roc_dir_name = renderText(paste(directory))
  
  baseline = NULL
  loaded_data = NULL
  active = NULL
  
  observeEvent(input$rocker_dir, {
    tryCatch({
      directory <<- choose_directory()
    },
    error = function(cond){
      directory <<- "No directory selected. Try again?"
      return(directory)
    })
    
    tryCatch({
      setwd(directory)
    },
    error = function(cond){
      directory <<- "No directory selected. Try again?"
      return(directory)
    })
    
    if(length(directory) == 0){
      directory <<- "No directory selected. Try again?"
    }
    
    
    if(directory == "No directory selected. Try again?"){
      output$roc_dir_name = renderText(paste(directory))
    }else{
      output$roc_dir_name = renderText(paste("Working in:", directory))
    }
    
  })
  
  observeEvent(input$load_dir, {
    tryCatch({
      descr = initial_data_gathering(directory)
      loaded = load_reads_and_describe(descr)
      
      baseline <<- descr
      loaded_data <<- loaded
      
    },
    error = function(cond){
      shinyalert("Error!", "Couldn't load ROCker Directory")
    })
    
    valid_prots = unique(loaded$origin)
    cats = descr$category[match(valid_prots, descr$protein_name)]
    updateCheckboxGroupInput(session, "proteins", label="Proteins to Show", choiceNames = paste(cats, ":", valid_prots) , choiceValues = valid_prots, selected = valid_prots)
    
  })
  
  observeEvent(input$proteins, {
    
    #print(input$proteins)
    
    if(!is.null(loaded_data)){
        active = loaded_data[origin %in% input$proteins ,]
        
        output$rocker_plot <- renderPlotly({
          
          if(!is.null(active)){
            colorlabs = c(Target = "blue3", Non_Target = "lightblue", Negative = "lightcoral")
            
            
            plot = ggplot(active, aes(x = sstart, xend = send, y = bitscore, yend = bitscore, color = classifier, label = origin_genome))+
              geom_segment()+
              xlim(c(0, upper))+
              ylim(c(0, top_bs))+
              xlab("Position in Protein")+
              ylab("Bitscore") +
              scale_color_manual(values = colorlabs)
            
            plot = ggplotly(plot)  
          }else{
            plot = warning_plot()
            plot = ggplotly(plot)
          }
          
          return(plot)
          
        })  
     }
    
  })
  

  
  session$onSessionEnded(function() {
    
    stopApp()
    
  })
  
}

runApp(list(ui = rocker_ui(), server = rocker_server), launch.browser = T)
