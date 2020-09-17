if (TRUE) {    # header
  #/*soh*************************************************************************
  # CODE NAME             : server.r
  # CODE TYPE  						: Program 
  # DATE OF UPDATE:         1-Feb-2019
  # DESCRIPTION           : Server code for BEACH app 
  # SOFTWARE/VERSION#     : R 3.3.0
  # INFRASTRUCTURE        : MS WINDOWS XP
  # LIMITED-USE MODULES   : N/A
  # BROAD-USE MODULES     : N/A
  # INPUT                 : 													
  # OUTPUT                : 
  #   VALIDATION LEVEL      : 
  #   REPLICATION PROGRAM   : N/A
  #   REGULATORY STATUS     : N/A
  #   REQUIREMENTS LOCATION : N/A
  #  -----------------------------------------------------------------------------
  #  Requirements: 
  #  -----------------------------------------------------------------------------
  #  Ver#  Author &                   Program History Description
  #        Peer Reviewer
  #  ----  ---------------            --------------------------------------------
  #  001   Danni Yu                   program (2013-now)
  #  002   Chenchen Yu                program (2013-2014)
  #  -----------------------------------------------------------------------------
}
if(!file.exists(local.path3)){dir.create(local.path3)}
myVersionCtr <<- "BEACH1.3.3" #2019-2-1 #add pump up window for plot

# 
# library(shinyAce)
# library(shiny)
# library(shinyjs)
# library(shinyBS)
# library(shinythemes)
# library(shinyWidgets)
# library(shinyFiles)
# library(shinyAce)
# library(shinycustomloader)
# 
# # htmlwidgets
# library(htmltools)
# library(htmlwidgets)
# library(formattable)
# library(kableExtra)
# library(methods); ##for get_data_demoUI
# # 
# # # viz and report packages
#  library(DT)
#  library(plotly)
#  library(visdat)
#  library(datadigest)
#  library(glue)


BeachUI <- fluidPage(
  
  #------- CSS code: Formats of radiobutton & checkbox -------#
  headerPanel(windowTitle='BEACH',
              tags$head(
                tags$link(rel = "stylesheet", type = "text/css"), #href = "bootstrap_dy.css"),
                tags$style(type="text/css","label.radio { display: inline-block; }",
                           "input[type=\'file\']{color: transparent;}",
                           ".radio input[type=\"radio\"] { float: none; }",
                           "label.checkbox { display: inline-block; }",
                           ".checkbox input[type=\"checkbox\"] { float: none; }",
                           #".shiny-output-error { visibility: hidden; }", #suppress red errors in shiny, to enable using visible
                           #".shiny-output-error:before { visibility: hidden; }",
                           "#TFL{max-width:95%; overflow-x:scroll; max-height:800px; overflow-y:scroll;}",
                           "#Input_outExpert{max-width:95%; overflow-x:scroll; max-height:400px; overflow-y:scroll;}",
                           "#widgetSide{max-width:100%; overflow-x:scroll;  max-height:100%; overflow-y:scroll;}",
                           "#add_analysis{height: 36px;}",
                           "#rcode{max-height:1000px; overflow-y:scroll;}",
                           "#scode{max-height:1000px; overflow-y:scroll;}",
                           "#pumpPlotOut{overflow-x:scroll; overflow-y:scroll;}",
                           ".shiny-text-output{max-width:95%; max-height:40%; overflow-y:scroll; }",
                           "th{text-align: center; border: 1px solid black;}",
                           "td{text-align: center; }", 
                           "table{border: 2px solid black;}"
                ))),
  
  #-----Start of Main Panel-----#
  absolutePanel( #for Analysis panel
    style = "background-color: #FFFFFF;",
    id = "controls1", class = "panel panel-default", fixed = FALSE,
    draggable=FALSE, top = "0%", left = "1%", right = "20%", bottom ="auto",
    width = "auto", height = "auto" , 
    uiOutput('tabs'),
    uiOutput('AnalysisTab'),
    uiOutput('SpecialTab') 
  ),
  
  #-----Start of Widget Panel-----#
  absolutePanel( #for widget/parameter panel
    style = "background-color: #F4F4F4;z-index: 200;",
    id = "controls1", class = "panel panel-default", fixed = FALSE,
    draggable = TRUE, top = "0%", left = "auto", right = "0%", bottom ="auto",
    width = "20%", height = "auto" , 
    p(style="color:black;",strong(myVersionCtr)),
    checkboxInput("collSidebar", "Show the sidebar for data input", value=TRUE),
    checkboxInput("useDT", "use renderDataTable", value=TRUE),
    checkboxInput("pumpPlot", "Show plot in a single window", value=FALSE),
    # checkboxInput("rundPAI", "run dPAI", value=FALSE),
    radioButtons("ncol.widg.rd", "number of columns", c(2, 1, 3), inline=TRUE),
    sliderInput('wpW', 'Width of the widget panel', 
                min=20, max=100, value=30, animate=TRUE),
    uiOutput('beachColor')
  ),  
  uiOutput('wp.width'),
  
  #-----Start of Sidebar Panel for LOA input-----#  
  absolutePanel( #for LOA panel
    style = "background-color: #F4F4F4;z-index: 500;",
    id = "controls2", class = "panel panel-default", fixed = FALSE,
    draggable = TRUE, top = "0%", left = "100%", right = "auto", bottom ="auto",
    width = "25%", height = "auto" , 
    
    
    #------LOA------# 
    p(style="color:black;",strong("List of Analysis")),
    fileInput('upfile', label='upload a LOA csv file',
              accept=c('.csv', 'text/csv', 
                       'text/comma-separated-values,text/plain'), multiple=FALSE),
    actionButton("load_analysis","load LOA"),
    downloadButton("save_loa","save LOA"),
    downloadButton("save_loa_cd", "save_CD"),
    checkboxInput('delete_loa', 'Only Selected', value=FALSE),
    wellPanel( uiOutput('LOA'), width="20%" ),
    radioButtons('landscp', 'RTF page layout', 
                 choices=c("Landscape", "Portrait"), selected="Portrait", inline=TRUE),
    radioButtons('onefileRTF', 'RTF output', 
                 choices=c('one file','multiple files'), selected='one file', inline=TRUE),
    downloadButton("save_output","OutputRTF"), #output as a RTF file    
    br(),
    textInput('hPath', label='HTM file location:', value=htmlPath), 
    actionButton("outputH", "Overwrite the HTM file"),
    uiOutput('save_outputH'),
    
    checkboxInput('expert', 'expert', value=FALSE),
    uiOutput('userExp'),
    checkboxInput('usage', 'usage', value=FALSE),
    checkboxInput('showAllSource', 'Show all R source code', value=FALSE),
    
    
    #----Data Manipulation, No UI output-----#
    uiOutput("getData"),
    uiOutput("getWidgets"),
    uiOutput("loa"),
    uiOutput('DTsel'),
    uiOutput('getDataSelected'), #defined to data_selected
    uiOutput("getTFL")
  ),    
  
  #-----Start of Sidebar Panel for data input-----#  
  div('data-display-if'="input.collSidebar == true", 
      absolutePanel( #for data input panel
        style = "z-index: 999;",
        id = "controls1", class = "panel panel-default", 
        fixed = FALSE,
        draggable = TRUE, top = "20%", left = "auto", right = "0%", bottom ="auto",
        width = "20%", height = "auto" , 
        
        #------Logo------#
        wellPanel(style = "background-color: #dbdbdb;", #D52B1E Lilly Red 100%
                  uiOutput("setTitle"),
#helpText(a('feedback', 
#target="_blank",
#href="https://gist.github.com/DanniYuGithub/74d2b327dcd63ff8e91633340ff95afa#file-forum-for-beach-users")),
                  #------Upload dataset------#
                  checkboxInput('data_reload', 'Make csv files reloadable with risk of messing up data.'),
                  uiOutput('status'), 
                  selectInput('comm_chr', "Comment Charactor in a CSV file", c("", "#", "*", "!", "~", "%", "^", "&")),
                  uiOutput('file'),
                  uiOutput("setconfig"),
                  selectInput('config.sel', 'Select A CD in the pool', choices=cdpool2),      
                  uiOutput('subsetcode1'), actionButton("submitcode","Submit Code"),
                  div(class='row'),
                  uiOutput('subsetcode2'),
                  actionButton("delete_tmpfile",'Clean trash on server.'), uiOutput('delete')
                  ))
      ),
  #-----End of Sidebar Panel-----#
  
  
  #-----Start pumpPlot single Window-----#  
  div('data-display-if'="input.pumpPlot == true", 
      absolutePanel( #for pumping out a single window
        style = "z-index: 999; background-color: #dbdbdb;",
        id = "controls1", class = "panel panel-default",  fixed=FALSE,
        draggable=TRUE, top="1%", left="20%", right="40%", bottom ="auto",
        width = "50%", height = "auto" , 
        plotOutput('pumpPlotOut')
        ))
  #-----End of Sidebar Panel-----#
  
  
  
) #Generate Tabs

shinyUI(BeachUI)