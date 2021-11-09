if(!require(zoo)){install.packages("zoo"); library(zoo)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(rjson)){install.packages("rjson"); library(rjson)}

# Download file from api

carregar_dados_drs <- function(estado = "SP", drs = "Grande_Sao_Paulo", data.base,
                               diff = FALSE, return = "list") {
  
  if (missing(data.base)){
    url = paste0("https://api.github.com/repos/covid19br/covid19br.github.io/contents/dados/DRS/", estado, "/", drs, "/tabelas_nowcasting_para_grafico")
    json_data <- fromJSON(paste(readLines(url, warn = FALSE), collapse=""))
    json_data_frame <- as.data.frame(json_data)
    # Filter objects
    name_columns <- grepl("name", colnames(json_data_frame))
    out <- json_data_frame[,name_columns]
    
    now_columns <- grepl("nowcasting_diario_srag_20", unlist(out[1,]))
    out2 <- out[,now_columns]
    
    # Extract last date
    datal <- gsub(".csv.*", "", gsub(".*srag_", "", unlist(out2[1,])))
    data.base <- max(datal)
  }
  
  nome.dir <- paste0("https://github.com/covid19br/covid19br.github.io/raw/master/dados/DRS/", estado, "/", drs, "/tabelas_nowcasting_para_grafico/")
  
  #data.base <- "2020_08_21"
  
  ######## CASOS ######## 
  ### covid ###
  now.covid <- read.csv(paste0(nome.dir,"nowcasting_diario_covid_",data.base,".csv"))
  now.covid.zoo <- cumsum(zoo(now.covid$estimate.merged, as.Date(now.covid$data)))
  # corta início e últimos 5 dias
  now.covid.zoo <- window(now.covid.zoo, start=as.Date('2020-3-1'))
  
  ### SRAG ###
  now.srag <- read.csv(paste0(nome.dir,"nowcasting_diario_srag_",data.base,".csv"))
  now.srag.zoo <- cumsum(zoo(now.srag$estimate.merged, as.Date(now.srag$data)))
  # corta início e últimos 5 dias
  now.srag.zoo <- window(now.srag.zoo, start=as.Date('2020-3-16'))
  
  ######## ÓBITOS ######## 
  
  ### covid ###
  now.obito.covid <- read.csv(paste0(nome.dir,"nowcasting_diario_obitos_covid_",data.base,".csv"))
  now.obito.covid.zoo <- cumsum(zoo(now.obito.covid$estimate.merged, as.Date(now.obito.covid$data)))
  # corta início e últimos 5 dias
  now.obito.covid.zoo <- window(now.obito.covid.zoo, start=as.Date('2020-3-1'))
  
  ### SRAG ###
  now.obito.srag <- read.csv(paste0(nome.dir,"nowcasting_diario_obitos_srag_",data.base,".csv"))
  now.obito.srag.zoo <- cumsum(zoo(now.obito.srag$estimate.merged, as.Date(now.obito.srag$data)))
  # corta início e últimos 5 dias
  now.obito.srag.zoo <- window(now.obito.srag.zoo, start=as.Date('2020-3-16'))
  
  if(return == "zoo") { 
    zoo.list <- merge.zoo(now.srag.zoo = now.srag.zoo,
                          now.obito.srag.zoo = now.obito.srag.zoo,
                          now.covid.zoo = now.covid.zoo,
                          now.obito.covid.zoo = now.obito.covid.zoo)
    if(diff == TRUE) zoo.list <- zoo.list %>% apply(2, diff)
  } 
  
  if(return == "list") { 
    
    zoo.list <- list(now.srag.zoo = now.srag.zoo,
                     now.obito.srag.zoo = now.obito.srag.zoo,
                     now.covid.zoo = now.covid.zoo,
                     now.obito.covid.zoo = now.obito.covid.zoo,
                     data.base = data.base)
    return(zoo.list)
  }
  
  return(zoo.list)
  
}
