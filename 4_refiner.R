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

protein_names = unlist(strsplit(protein_names, split = "_read_length_"))

read_lens = protein_names[c(F,T)]
protein_names = protein_names[c(T,F)]

df[, protein_name := protein_names]
df[, mean_length := as.numeric(read_lens)]
rm(protein_names)

return(df)
}


purity = function(vector){
  counts = table(vector)
  
  return(max(counts)/sum(counts))
  
}

collect_ma = function(directory){
  
  files = list.files(full.names = T, path = directory, pattern = "_MA.txt", recursive = T)
  if(length(files)== 0){
    return(NULL)
  }else{
    ma = files[grepl("_MA.txt", files)]
    ma = readLines(ma)
    
    defline_lines = grepl(">", ma)
    deflines = ma[defline_lines]
    #Clean out deflines
    ma = ma[!defline_lines]
    
    #Find how many rows are sequence and what they go with.
    grps = rle(defline_lines)$lengths[c(F,T)]
    
    #Collect lines by sequence
    grab_these = rep(seq(1:length(grps)), times = grps)
    
    grouped_seqlines = split(ma, f = grab_these)
    grouped_seqlines = unname(grouped_seqlines)
    grouped_seqlines = lapply(grouped_seqlines, function(x){
      joined = paste(x, collapse = '')
      splits = strsplit(joined, split = '')[[1]]
    })
    
    prot_ids = unlist(lapply(deflines, function(x){
      name = strsplit(x, split = "\t")[[1]][1]
      name = substr(name, 2, nchar(name))
      }))

    df = data.table(defline = rep(deflines, times = lengths(grouped_seqlines)), prot_id = rep(prot_ids, times = lengths(grouped_seqlines)), position = rep(1:length(grouped_seqlines[[1]]), times = length(grouped_seqlines)), base = unlist(grouped_seqlines))

    return(df)
    
  }
  
}

collect_ma_seqs <- function(directory){
  files = list.files(full.names = T, path = directory, pattern = "_MA.txt", recursive = T)
  if(length(files)== 0){
    return(NULL)
  }else{
    ma = files[grepl("_MA.txt", files)]
    ma = readLines(ma)
    
    defline_lines = grepl(">", ma)
    deflines = ma[defline_lines]
    #Clean out deflines
    ma = ma[!defline_lines]
    
    #Find how many rows are sequence and what they go with.
    grps = rle(defline_lines)$lengths[c(F,T)]
    
    #Collect lines by sequence
    grab_these = rep(seq(1:length(grps)), times = grps)
    
    grouped_seqlines = split(ma, f = grab_these)
    grouped_seqlines = unname(grouped_seqlines)
    grouped_seqlines = unlist(lapply(grouped_seqlines, function(x){
      joined = paste(x, collapse = '')
    }))
    
    prot_ids = unlist(lapply(deflines, function(x){
      name = strsplit(x, split = "\t")[[1]][1]
      name = substr(name, 2, nchar(name))
    }))
    
    df = data.table(deflines, prot_ids, grouped_seqlines)
    
    return(df)
    
  }
  
}


#Use the file paths to load the files
load_reads_and_describe = function(rdf){
  
    read_data = rbindlist(mapply(function(x, y){
      
      protein_name = unlist(strsplit(x, split = "/aligned_reads/"))[c(F,T)]
      protein_name = unlist(strsplit(protein_name, split = "_aligned_reads.txt"))
      
      protein_name = unlist(strsplit(protein_name, split = "_read_length_"))
      
      #read_lens = protein_name[c(F,T)]
      protein_name = protein_name[c(T,F)]
      
      
      if(file.size(x) == 0){
        #Suppress error message
        one_file = data.table(NULL)
      }else{
        
        one_file = fread(x, header = F, )
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
          one_file[, read_length := y]
          
        }
      }
      
      return(one_file)
      
    }, rdf$path, rdf$mean_length))
 
    read_data$read_length = factor(read_data$read_length, levels = sort(unique(read_data$read_length)))
    
    read_data = split(read_data, read_data$read_length)
  
  return(read_data)
  
}

#rows = input$vertical_resolution
#sliding_window = input$sliding_window
collect_roc_curve = function(mapped_reads_df, rows = 50, sliding_window = 20, upper, top_bs, low_bs){
  
  #If the summarization checkbox isn't clicked, we bin the reads and show that plot instead.
  
  #rows = input$vertical_resolution
  
  #Set up summary datatable
  summary_dt = data.table(CJ(seq(0, upper, 1), 0:rows), on_target = 0, confounder = 0, ratio = 0)
  colnames(summary_dt) = c("x", "y", "on_target", "confounder", "ratio")
  #bin the results
  step = max(mapped_reads_df$bitscore)/rows
  
  #figure out the y bin for each read
  mapped_reads_df[, ybin := bitscore %/% step]
  
  
  split_vec = mapped_reads_df$classifier == "Target"
  #If either is empty, we get a null data table, which is a case handled by the below ifs
  #Select hits
  tgt = mapped_reads_df[split_vec,]
  #Select misses
  neg = mapped_reads_df[!split_vec]
  
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
  roc_data = data.table(wstart = 1:(upper - sliding_window), wend = (sliding_window+1):upper, most_discriminant = NA)
  
  for(i in 1:(upper - sliding_window)){
    sub = summary_dt[Pos_in_Protein >= i & Pos_in_Protein <= i + sliding_window, list(tgt = max(Target), conf = max(Confounder)), by = Bitscore]
    
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
  
  #roc_model <<- roc_data
  
  #Fill out the data so that the lines are pretty.
  first = roc_data[1,]
  first$midpt = 0
  last = roc_data[nrow(roc_data),]
  last$midpt = last$wend
  roc_data = rbindlist(list(first, roc_data, last))
  
  #Clean up
  mapped_reads_df[, ybin := NULL, ]
  
  return(list(summary_dt, roc_data))
  
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
  tabsetPanel(id = "tabs",
  
  tabPanel("ROCkOut Model",              
  sidebarPanel(width = 3,
         #Selectable options go here
         
         actionButton('rocker_dir', 'Select a ROCker project', icon = icon("folder-open")),
         verbatimTextOutput("roc_dir_name"),
         
         actionButton('load_dir', "Load the ROCker directory", icon = icon("folder-upload")),
         
         sliderInput('sliding_window', 'Rolling Average', value = 20, max = 150, min = 0),
         
         sliderInput('vertical_resolution', 'ROC Resolution', value = 50, max = 100, min = 0),
         
         selectInput('select_read_length', 'Read Length Selection', choices = c("Load a project first" = "load"), selected = NULL),
         
         actionButton('print_model', 'Output ROCkOut Model', icon = icon('dollar-sign')),
         
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
  
  #model_vis_tab end
  ),
  
  tabPanel("3D Plot",
          sidebarPanel(width = 3,
                       checkboxInput("roc_surface", "Show ROCkOut Surface?", value = FALSE)
                       ),
           mainPanel(id = "rotation_plot",
                     fluidRow(
                       column(12, align="center",
                       plotlyOutput("plot_3d", height = "850px")
                       )
                     )
             
           #End main panel
           )
           #End 3d tab
           ),
  
  tabPanel("Multiple Alignment",
           
           sidebarPanel(width=3,
                        checkboxInput("conservation", "Color by pct. conservation?", value = FALSE)),
           
           mainPanel(id = "multiple_aln",
                     #Spacing
                     fluidRow(
                       column(12,
                              div(style = "height:60px;background-color: white;", "")
                       )
                     ),
                     fluidRow(
                       
                       plotOutput("ma_plot", height = "850px")
                       
                     )
                     )
           
           )
  
  #tabset end parentheses
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
  roc_model = NULL
  
  target_seqs = NULL
  active_seqs = NULL
  
  ma_data = NULL
  active_ma = NULL
  
  tree_data = NULL
  active_tree = NULL
  
  plot_3d = NULL
  
  #For use in the plotting functions of the server.
  colorlabs = c(Target = "blue3", Non_Target = "lightblue", Negative = "lightcoral")
  
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
      
      tryCatch({
        mult_aln = collect_ma(directory)
        ma_data <<- mult_aln
        
      }, error = function(cond){
        shinyalert("Error!", "Multiple alignment data not found!\n\nROCkOut will continue, but you will not be able to use the multiple alignment tab.")
      })
      
      baseline <<- descr
      loaded_data <<- loaded

      valid_prots = unique(unlist(lapply(loaded, function(x){
        return(unique(x$origin))
      })))
      
      valid_lengths = as.character(sort(unique(descr$mean_length)))
      names(valid_lengths) = paste0("Average Read Length: ", valid_lengths)
      
      cats = descr$category[match(valid_prots, descr$protein_name)]
      
      
      updateCheckboxGroupInput(session, "proteins", label="Proteins to Show", choiceNames = paste(cats, ":", valid_prots) , choiceValues = valid_prots, selected = valid_prots)
      
      
      updateSelectInput(session, "select_read_length", choices = valid_lengths, selected = valid_lengths[1])
      
    },
    error = function(cond){
      shinyalert("Error!", "Couldn't load ROCker Directory")
    })
    

    
  })
  
  observeEvent(input$proteins, {
    
    if(!is.null(loaded_data)){
      
      index = which(names(loaded_data) == input$select_read_length)
      active <<- loaded_data[[index]]
      
      active <<- active[origin %in% input$proteins, ]
      
        output$rocker_plot <- renderPlotly({
          
          if(!is.null(active)){
            
            upper = max(active$send)
            top_bs = max(active$bitscore)
            low_bs = min(active$bitscore)
            model_data = collect_roc_curve(active, input$vertical_resolution, input$sliding_window, upper, top_bs, low_bs)
            
            summary_dt = model_data[[1]]
            roc_data = model_data[[2]]
            
            model_data = NULL
            
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
    
    if(!is.null(ma_data)){
     active_ma <<- ma_data[prot_id %in% input$proteins,] 
     
     output$ma_plot <- renderPlot({
     
     if(is.null(active_ma)){
       ma_plot = warning_plot()
     }else{
       
       if(input$conservation){
         active_ma[, conservation := purity(base), by = position]
         ma_plot = ggplot(active_ma, aes(x = position, y = prot_id, fill = conservation))+
           geom_tile() + 
           ylab("Protein Sequence")+
           xlab("AA in Protein")
         
       }else{
         
        #Make sure non-alns are blank
        if("-" %in% unique(active_ma$base)){
          group.colors = c("grey65", scales::hue_pal()(length(unique(active_ma$base))-1))
          others = sort(unique(active_ma$base))
          others = others[!others %in% '-']
          names(group.colors) = c("-", others)
        }else{
          #this has no non-aln
          group.colors = c(scales::hue_pal()(length(unique(active_ma$base))))
          others = sort(unique(active_ma$base))
          names(group.colors) = c(others)
        }

         ma_plot = ggplot(active_ma, aes(x = position, y = prot_id, fill = base))+
           geom_tile() + 
           ylab("Protein Sequence")+
           xlab("AA in Protein") +
           scale_fill_manual(values = group.colors)
       }
     }
     
     return(ma_plot)
     
     })
     
     
    }
    
  })
  
  observeEvent(input$select_read_length, {
    if(!is.null(loaded_data)){
      
      index = which(names(loaded_data) == input$select_read_length)
      active <<- loaded_data[[index]]
      
      output$rocker_plot <- renderPlotly({
        
        if(!is.null(active)){
          
          upper = max(active$send)
          top_bs = max(active$bitscore)
          low_bs = min(active$bitscore)
          model_data = collect_roc_curve(active, input$vertical_resolution, input$sliding_window, upper, top_bs, low_bs)
          
          summary_dt = model_data[[1]]
          roc_data = model_data[[2]]
          
          model_data = NULL
          
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
  
  
  observeEvent(input$tabs, {
    
    if(input$tabs == "3D Plot"){
      output$plot_3d <- renderPlotly({
        
        if(!is.null(active)){
          
          three_dee = rbindlist(loaded_data)
          three_dee = three_dee[origin %in% input$proteins, ]
          
          upper = max(three_dee$send)
          top_bs = max(three_dee$bitscore)
          low_bs = min(three_dee$bitscore)
          
          three_dee[, read_length := as.numeric(as.character(read_length)),]
          
          
          plot_3d <<- plot_ly(three_dee, x = ~abs((sstart+send)/2), y = ~read_length, z = ~bitscore, color = ~classifier) %>% add_markers()
          
          
        }else{
          plot = warning_plot()
          plot = ggplotly(plot)
        }
        
        return(plot_3d)
        
      })  
    }
    
  })
  
  observeEvent(input$roc_surface, {
    
    if(!is.null(plot_3d)){
      
      roc_curves = lapply(loaded_data, function(x){
        
        upper = max(x$send)
        top_bs = max(x$bitscore)
        low_bs = min(x$bitscore)
        
        rd = collect_roc_curve(x, input$vertical_resolution, input$sliding_window, upper = upper, top_bs = top_bs, low_bs = low_bs)
        
        just_curve = rd[[2]]
        
        just_curve$read_length <- as.numeric(as.character(x$read_length[1]))
        
        return(just_curve)
        
      })
      
      combined_curves = rbindlist(roc_curves)
      
      grouped = split(combined_curves, f = combined_curves$midpt)
      
      rm(combined_curves)
      
      surface = rbindlist(mapply(function(x, mid){
        model = lm(x$most_discriminant ~ poly(x$read_length, 2, raw = TRUE))
        model = model$coefficients
        
        space = seq(min(x$read_length), max(x$read_length), by = 1)
        
        preds = space * model[2] + (space^2 * model[3]) + model[1]
        
        slice = data.table(x = mid, y = space, z = preds)
        
        #print(slice)
        
        return(list(slice))
        
      }, grouped, as.numeric(names(grouped))))
      
      
      
      
    }
    
  })
  
  #This will need an update.
  observeEvent(input$print_model, {
    
    if(!is.null(roc_model)){
      dirname = paste0(directory, "/ROCkOut_Model")
      
      #print(dirname)
      
      if(!dir.exists(dirname)){
        dir.create(dirname)
      }
    }
    
    
    model_output = paste0(dirname, "/ROCkOut_Model.txt")
    
    fwrite(roc_model, model_output, sep = "\t")
    
    if(!is.null(active_ma)){
      
      joint = collect_ma_seqs(directory)
      joint = joint[prot_ids %in% active_ma$prot_id,]
      
      for(i in 1:nrow(joint)){
        cat(joint$deflines[i], file=model_output, append=TRUE, sep = "\n")
        cat(joint$grouped_seqlines[i], file=model_output, append=TRUE, sep = "\n")
      }
    }
    
  })
  
  session$onSessionEnded(function() {
    
    stopApp()
    
  })
  
}

runApp(list(ui = rocker_ui(), server = rocker_server), launch.browser = T)
 