if(!require(zoo)){install.packages("zoo"); library(zoo)}
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)}
if(!require(rjson)){install.packages("rjson"); library(rjson)}

# Download file from api
carregar_dados_estado <- function(estado = "SP", data.base, diff = FALSE,
                                  return = "list", type.data = "nowcasting") {

    if(type.data == "nowcasting") {
        type.data <- "nowcasting_diario"
    } else {
        type.data <- "r_efetivo"
    }

    if (missing(data.base)) {
        url = paste0("https://api.github.com/repos/covid19br/covid19br.github.io/contents/dados/estado/", estado, "/tabelas_nowcasting_para_grafico")
        json_data <- fromJSON(paste(readLines(url, warn = FALSE), collapse=""))
        json_data_frame <- as.data.frame(json_data)
        # Filter objects
        name_columns <- grepl("name", colnames(json_data_frame))
        out <- json_data_frame[,name_columns]
        
       
        now_columns <- grepl(paste0(type.data, "_srag_20"), unlist(out[1,]))
        out2 <- out[,now_columns]
        
        # Extract last date
        datal <- gsub(".csv.*", "", gsub(".*srag_", "", unlist(out2[1,])))
        data.base <- max(datal)
    }
    
    nome.dir <- paste0("https://github.com/covid19br/covid19br.github.io/raw/master/dados/estado/", estado, "/tabelas_nowcasting_para_grafico/")
    
    #data.base <- "2020_08_21"
    
    ######## CASOS ######## 
    ### covid ###
    now.covid <- read.csv(paste0(nome.dir,type.data,"_covid_",data.base,".csv"))
    if (type.data == "nowcasting_diario"){
        now.covid.zoo <- cumsum(zoo(now.covid$estimate.merged, as.Date(now.covid$data)))
    } else {
        now.covid.zoo <- zoo(now.covid$Mean.R, as.Date(now.covid$data))
    }
    # corta início e últimos 5 dias
    now.covid.zoo <- window(now.covid.zoo, start=as.Date('2020-3-15'))
    
    ### SRAG ###
    now.srag <- read.csv(paste0(nome.dir,type.data,"_srag_",data.base,".csv"))
    if (type.data == "nowcasting_diario"){
        now.srag.zoo <- cumsum(zoo(now.srag$estimate.merged, as.Date(now.srag$data)))
    } else {
        now.srag.zoo <- zoo(now.srag$Mean.R, as.Date(now.srag$data))
    }
    # corta início e últimos 5 dias
    now.srag.zoo <- window(now.srag.zoo, start=as.Date('2020-3-15'))
    
    ######## ÓBITOS ######## 
    if (type.data == "nowcasting_diario"){
        ### covid ###
        now.obito.covid <- read.csv(paste0(nome.dir,type.data,"_obitos_covid_",data.base,".csv"))
        now.obito.covid.zoo <- cumsum(zoo(now.obito.covid$estimate.merged, as.Date(now.obito.covid$data)))
        # corta início e últimos 5 dias
        # now.obito.covid.zoo <- window(now.obito.covid.zoo, start=as.Date('2020-3-15'))
        
        ### SRAG ###
        now.obito.srag <- read.csv(paste0(nome.dir,type.data,"_obitos_srag_",data.base,".csv"))
        now.obito.srag.zoo <- cumsum(zoo(now.obito.srag$estimate.merged, as.Date(now.obito.srag$data)))
        # corta início e últimos 5 dias
        # now.obito.srag.zoo <- window(now.obito.srag.zoo, start=as.Date('2020-3-15'))
    }
    
    
    if(return == "zoo") { 
        if (type.data == "nowcasting_diario"){
            zoo.list <- merge.zoo(now.srag.zoo = now.srag.zoo,
                                  now.obito.srag.zoo = now.obito.srag.zoo,
                                  now.covid.zoo = now.covid.zoo,
                                  now.obito.covid.zoo = now.obito.covid.zoo)
        } else {
            zoo.list <- merge.zoo(reff.srag.zoo = now.srag.zoo,
                                  reff.covid.zoo = now.covid.zoo)
            
        }
        if(diff == TRUE) zoo.list <- zoo.list %>% apply(2, diff)
        
    } 
    
    if(return == "list") { 
        if(type.data == "nowcasting_diario"){
            zoo.list <- list(now.srag.zoo = now.srag.zoo,
                             now.obito.srag.zoo = now.obito.srag.zoo,
                             now.covid.zoo = now.covid.zoo,
                             now.obito.covid.zoo = now.obito.covid.zoo)
        } else {
            zoo.list <- list(reff.srag = now.srag.zoo,
                             reff.covid = now.covid.zoo)
        }
        
        return(zoo.list)
    }
    
    return(zoo.list)
    
}
