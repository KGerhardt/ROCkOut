setwd("C:/Users/Kenji/Desktop/rocker_proj/roc_out")

library(ggplot2)
library(data.table)
library(plotly)
library(reticulate)
library(shiny)
library(shinyalert)
library(shinyBS)



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
    
    if(file.size(x) == 0){
      #Suppress error message
      one_file = data.table(NULL)
    }else{
    
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
         
         sliderInput('sliding_window', 'Rolling Average', value = 20, max = 150, min = 0),
         
         sliderInput('vertical_resolution', 'ROC Resolution', value = 50, max = 100, min = 0),
         
         checkboxInput('summarize', 'Summarize Reads', value = TRUE),
         bsTooltip('summarize', 'If this box is checked, reads will be displayed in bins according to the ROC Resolution slider. This keeps plots fast and responsive even with many reads.'),
         
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
    
    #Switch to summary if enough reads - 30k is a decent point.
    if(nrow(loaded_data) > 3e04){
      shinyalert(title = "Summarizing", text = paste0("There were ", nrow(loaded_data), " reads. Switching to summary plot."), type = 'info')
      updateCheckboxInput(session, "summarize", value = FALSE)
    }
    
  })
  
  observeEvent(input$proteins, {
    
    #print(input$proteins)
    
    if(!is.null(loaded_data)){
        active = loaded_data[origin %in% input$proteins ,]
        
        output$rocker_plot <- renderPlotly({
          
          if(!is.null(active)){
            colorlabs = c(Target = "blue3", Non_Target = "lightblue", Negative = "lightcoral")
            
            upper = max(active$send)
            top_bs = max(active$bitscore)
            low_bs = min(active$bitscore)
            
            #If the summarization checkbox isn't clicked, we bin the reads and show that plot instead.
            
            rows = input$vertical_resolution
            
            #Set up summary datatable
            summary_dt = data.table(CJ(seq(0, upper, 1), 0:rows), on_target = 0, confounder = 0, ratio = 0)
            colnames(summary_dt) = c("x", "y", "on_target", "confounder", "ratio")
            #bin the results
            step = max(active$bitscore)/rows
            
            #figure out the y bin for each read
            active[, ybin := bitscore %/% step]
            
            
            split_vec = active$classifier == "Target"
            #If either is empty, we get a null data table, which is a case handled by the below ifs
            #Select hits
            tgt = active[split_vec,]
            #Select misses
            neg = active[!split_vec]
            
            rm(split_vec)
            
            if(nrow(tgt) > 0){
              #This is a fairly complicated expression. Break it down:
              #Split each observation by y bin - effectively, consider reads only within a specific bitscore window and perform each of the following for each bin:
              
              #create a vector counting from sstart to ssend, inclusive (positions a read covers)
              #unlist those vectors so that all of the covered positions are listed.
              #sort the positions so that repeat covers of the same positions are adjacent in the overall list
              #get a run-length encoding on those results, effectively a count of occurrences per position
              
              tgt = tgt[, rle(sort(unlist(mapply(":", sstart, send)))), by = ybin]
              tgt[, label := "On-Target"]
              
              #match x and y bins; assign count of reads covering each position at each y bin to corresp. summary DT bin.
              summary_dt[tgt, on_target := lengths, on = c(x = "values", y = "ybin")]
              
              #clean up
              rm(tgt)
            }
            if(nrow(neg) > 0){
              #repeat the process for negatives
              neg = neg[, rle(sort(unlist(mapply(":", sstart, send)))), by = ybin]
              neg[, label := "Confounder"]
              
              summary_dt[neg, confounder := lengths, on = c(x = "values", y = "ybin")]
              
              rm(neg)
            }
            
            #summary_dt = rbindlist(list(tgt, neg))
            #colnames(summary_dt) = c("Bitscore", "Count", "Pos. in Protein", "Group")
            
            #We don't need non-data, so remove it
            summary_dt = summary_dt[on_target > 0 | confounder > 0,]
            
            #0 to 1 for positive > negative
            summary_dt[on_target >= confounder, ratio := on_target/(confounder + on_target)]
            #-1 to 0 for negative > positive
            summary_dt[on_target < confounder, ratio := -confounder/(confounder + on_target)]
            
            #Convert back to the real bin values for y
            summary_dt[, y := y * step]
            summary_dt[, bpct := on_target + confounder]
            
            colnames(summary_dt) = c("Pos_in_Protein", "Bitscore", 'Target', 'Confounder', "Ratio", "BP_Count")
            
            #Caculate ROC overlay
            #Prepare a dataframe with sliding windows describing starts and ends along the whole protein.
            roc_data = data.table(wstart = 1:(upper - input$sliding_window), wend = (input$sliding_window+1):upper, most_discriminant = NA)
            

            for(i in 1:(upper - input$sliding_window)){
              sub = summary_dt[Pos_in_Protein >= i & Pos_in_Protein <= i + input$sliding_window, list(tgt = max(Target), conf = max(Confounder)), by = Bitscore]

              #Order descending
              setorder(sub, -Bitscore)
              
              #Cumulative sums...
              sub[, tgt := cumsum(tgt)]
              sub[, conf := cumsum(conf)]
              #We care about maximizing the Youden index, per the original ROCker paper
              #Youden = [TP / (TP + FN)] + [TN / (FP + TN)] - 1
              
              #TP = tgt
              #FN = max(tgt) - tgt
              #TN = max(conf) - conf
              #FP = conf
              
              #[TP / (TP + FN) = tgt / (tgt + max(tgt) - tgt) = [tgt / max(tgt)]
              #[TN / (FP + TN)] = max(conf) - conf / conf + max(conf) - conf = [max(conf) - conf / max(conf)]
              
              
              sub[, Youden := ((tgt / max(tgt)) + (max(conf) - conf / max(conf))) - 1, ]
              
              if(all(is.nan(sub$Youden))){
                cutoff = min(sub$Bitscore)
              }else{
                cutoff = sub[Youden == max(Youden, na.rm = T), Bitscore]
              }
              
              
              roc_data$most_discriminant[i] = cutoff
            }
            
            
            roc_data[, midpt:=(wstart + wend)/2]
            #Fill out the data so that the lines are pretty.
            first = roc_data[1,]
            first$midpt = 0
            last = roc_data[nrow(roc_data),]
            last$midpt = last$wend
            roc_data = rbindlist(list(first, roc_data, last))
            
            #####################
            
            if(input$summarize){
              #Use the summary data instead
              plot = ggplot(summary_dt, aes(x = Pos_in_Protein, y = Bitscore, fill = Ratio, label = BP_Count)) +
                geom_raster()+
                xlim(c(-1, upper+1))+
                ylim(c(low_bs-1, top_bs+1))+
                xlab("Position in Protein")+
                ylab("Bitscore") +
                scale_fill_gradient2(low = "lightcoral", mid = "grey65", high = "blue3", midpoint = 0, limits = c(-1, 1)) + 
                geom_step(data = roc_data, aes(x = roc_data$midpt, y =roc_data$most_discriminant), lwd = 1.5, inherit.aes = F)
              
            }else{
              
              #active is the per-read dataframe. This can be computationally difficult to display, so we toggle it off by default.
              plot = ggplot(active, aes(x = sstart, xend = send, y = bitscore, yend = bitscore, color = classifier, label = origin_genome))+
                geom_segment()+
                xlim(c(-1, upper+1))+
                ylim(c(low_bs-1, top_bs+1))+
                xlab("Position in Protein")+
                ylab("Bitscore") +
                scale_color_manual(values = colorlabs) + 
                geom_step(data = roc_data, aes(x = roc_data$midpt, y =roc_data$most_discriminant), lwd = 1.5, inherit.aes = F)
              
            }
              

            plot = ggplotly(plot)
              
            
            
            #Add in ROC overlay.
            
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
 