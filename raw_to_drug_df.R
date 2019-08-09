###
# By Oskar Gauffin and modified by Marcus Westerberg
# Modified: 2019-08-02

# !diagnostics off
# The line above prevents dplyr from outputting a lot of nonsense warnings

##############################################
# Line calculation and drug duration-script
##############################################
# This script is written to calculate crpc-line of treatments, and from this, treatment duration, 
# based on the cessation (uts?ttning) of the drug or death date of the patient.
# A few exceptions are built in.
# Docetaxel started within nine month of pre-GnRH-metastases (either from NPCR (at dx) or PPC (later)) in the first line, is not counted
# as a crpc-treatment.
# Other CT (CM, KEES, etc) could be multi-drug-line. Thus, secondary initiations of other CT only counts as a new line of treatment if some other 
# drug (not in "Other CT") is in between two initiations of "Other CT"-drugs. 

# This data set is not easy to clean. Initiation and cessation could occur at the same date, and then preceding cessations will decide
# whether this is a one-day-line, or an ended line and a new line started on the same date. 
# It was not considered sensible to require that the initiation is marked as an initiation, as some have free text comments or "side effects" as 
# initiation comment (probably, the initiation was due to some side effect of another drug, or possible entered retrospectively that the drug was initiated
# but with side effects found later.)
# Xofigo is sometimes mistakenly ended instead of paused after each cycle. We try to fix this. 
# The overall improvement compared to the previous definition (line was from first date to last date of the drug) is to find multiple lines for the same drug
# this was not the case in the previous script. That is, one day docetaxel, ten year abiraterone, one day docetaxel should be considered as three lines,
# not as a first line of docetaxel with duration 10 years and 2 days, and an overlapping "second" line of abiraterone.

### Convert raw-PPC to usable PPC-INCA-files.
### Data cleaning -----------------------------
nothing <- lapply(c("magrittr","lubridate","plyr","dplyr", "stringr","data.table"), library, character.only=T)

dst_base <- haven::read_sas(paste0("R:\\PCBase41\\Works41\\PPC\\2018\\_Xofigostudie\\Data_source\\PPC\\sos_ppc1_cleaned_dst_ba_lopnr.sas7bdat")) # base

dst_base %<>% select(LopNr, birth_dat, death_dat, LopNr, firstmet_dat) %>% rename(firstmet_dat_base = firstmet_dat)

## Add M1-data ----------------
# For PPC we only care if the metastasis occurred BEFORE initiation of GnRH
dst_met <- dst_base %>% select(LopNr, firstmet_dat_base) %>% filter(!is.na(firstmet_dat_base))


dst_diag <- haven::read_sas(paste0("R:\\PCBase41\\Works41\\PPC\\2018\\_Xofigostudie\\Data_source\\PPC\\sos_ppc3_cleaned_dst_di_lopnr.sas7bdat")) # diag

dst_diag$image_metast[dst_diag$image_metast %in% ""] <- NA
dst_diag$created_dat[dst_diag$created_dat %in% ""] <- NA
dst_diag %<>% filter(image_assess %in% c("Blandad Respons", "Progress") | !is.na(image_metast)) %>%
  group_by(LopNr) %>% mutate(first_met_ppc = min(image_dat)) %>% ungroup() %>% select(LopNr, first_met_ppc) %>% distinct()





# Find first initiation of GnRH
drugs <- haven::read_sas(paste0("R:\\PCBase41\\Works41\\PPC\\2018\\_Xofigostudie\\Data_source\\PPC_raw\\sos_drug_tab_lopnr.sas7bdat")) # diag #drug raw
#drugs[drugs==""] = NA
# drugs$LMIH_handelsedatum <- as.Date(drugs$LMIH_handelsedatum)
gnrh <- tolower(c("Leuprorelin","Goserelin","Buserelin","Triptorelin","Histrelin"))
drugs$LMIH_handelsedatum[drugs$LMIH_handelsedatum %in% "NA"] <- NA

gnrh_df <- drugs %>% rename(drug_name = sub_substans_vd_sub_subnamnrek_) %>% filter(drug_name %in% gnrh) %>% group_by(LopNr) %>% 
  mutate(first_gnrh = min(c(LMIH_handelsedatum,Inf), na.rm=T)) %>% select(LopNr, first_gnrh)

sum(gnrh_df$first_gnrh==Inf)
gnrh_df <- gnrh_df[!gnrh_df$first_gnrh==Inf,]



merge(dst_met, dst_diag, all.x=T) -> dst_met 
merge(dst_met, gnrh_df, all.x=T) -> dst_met 
base_met_after_gnrh <- !is.na(dst_met$first_gnrh) & !is.na(dst_met$firstmet_dat_base) & dst_met$first_gnrh <= dst_met$firstmet_dat_base
diag_met_after_gnrh <- !is.na(dst_met$first_gnrh) & !is.na(dst_met$firstmet_dat_ppc) & dst_met$first_gnrh <= dst_met$firstmet_dat_ppc
dst_met[base_met_after_gnrh, "firstmet_dat_base"] = NA
dst_met[diag_met_after_gnrh, "firstmet_dat_base"] = NA
dst_met %<>% select(-first_gnrh)




# For NPCR we use all metastases
npcr <- data.table::fread("R:/PCBase41/BasicData41/RCC/NPCR/csv_versioner/npcr.csv",
                          na.strings=c("", "NA", NA),
                          fill=T,
                          colClasses="character",
                          data.table=F,
                          showProgress=FALSE)

npcr %<>% select(LopNr, d0_date, d0_m_txt) %>% filter(d0_m_txt=="M1") %>% select(LopNr, first_met_npcr=d0_date, -d0_m_txt) %>% 
  mutate(first_met_npcr = substr(first_met_npcr, 0, 10))

merge(dst_met, npcr, all.x=T) -> dst_met
dst_met$first_met_npcr <- as.Date(dst_met$first_met_npcr)
dst_met$firstmet_dat_base <- as.Date(dst_met$firstmet_dat_base)

dst_met %<>% filter(any(!is.na(firstmet_dat_base)|!is.na(first_met_ppc)|!is.na(first_met_npcr)))

dst_met$firstmet_all = NA
dst_met %<>% distinct()
dst_met$firstmet_dat_base[dst_met$firstmet_dat_base %in% "NA"]<-NA
dst_met$first_met_ppc[dst_met$first_met_ppc %in% "NA"]<-NA
dst_met$first_met_npcr[dst_met$first_met_npcr %in% "NA"]<-NA

all_x <- split(dst_met, factor(dst_met$LopNr))

ldply(all_x, function(x){
  # x <- all_x[[5]] 
  # x <- all_x[[1]] 
  if(any(!is.na(x[,2:4]))){
    x$firstmet_all <- x[,2:4][which.min(as.numeric(x[,2:4]))][1,1]
  }
  x
}, .id="LopNr") -> dst_met 
dst_met %<>% select(LopNr, firstmet_all)

drugs <- haven::read_sas(paste0("R:\\PCBase41\\Works41\\PPC\\2018\\_Xofigostudie\\Data_source\\PPC_raw\\sos_drug_tab_lopnr.sas7bdat")) # diag #drug raw
#drugs[drugs==""] = NA
drugs <- merge(dst_base, drugs, all.y=T)
drugs <- drugs %>% rename(drug_name = sub_substans_vd_sub_subnamnrek_)
drugs$drug_name <- gsub(" \\(vattenfri\\)", "", drugs$drug_name)
#drugs <- drugs[drugs$birth_dat %in% "1934-02-26",]
drugs$LMIH_handelseorsak_Beskrivning[drugs$LMIH_handelseorsak_Beskrivning %in% "Nyins?ttning"] <- "Nyins?ttning"
drugs$LMIH_handelseorsak_Beskrivning[drugs$LMIH_handelseorsak_Beskrivning %in% "?terins?ttning"] <- "?terins?ttning"
drugs$LMIH_handelseorsak_Beskrivning[drugs$LMIH_handelseorsak_Beskrivning %in% "Patientens ?nskem?l"] <- "Patientens ?nskem?l"
drugs$LMIH_handelseorsak_Beskrivning[drugs$LMIH_handelseorsak_Beskrivning %in% "Str?lbehandling"] <- "Str?lbehandling"


start_var <- c("LopNr","drug_name", "LMIH_handelseorsak_Beskrivning","LMIH_handelse_Beskrivning",
               "LMIH_handelsedatum", "LMIH_handelseorsakann", "death_dat")
stop_var <- c("LopNr","drug_name", "LMIU_utsattningsdatum", "LMIU_utsattningsorsak_Beskrivni",
              "LMIU_utsattningsorsakann", "death_dat")

drug_df <- bind_rows(drugs[, start_var] %>% rename(event = LMIH_handelseorsak_Beskrivning,
                                                   comment=LMIH_handelse_Beskrivning,
                                                   event_dat=LMIH_handelsedatum,
                                                   free_text_comment=LMIH_handelseorsakann),
                     
                     drugs[, stop_var] %>% rename(comment=LMIU_utsattningsorsak_Beskrivni,
                                                  event_dat=LMIU_utsattningsdatum,
                                                  free_text_comment=LMIU_utsattningsorsakann)) %>% filter(!is.na(event_dat))

drug_df$event[is.na(drug_df$event)] = "Uts?ttning"
drug_df <- drug_df[!duplicated(drug_df),]
drug_df <- drug_df[order(drug_df$event_dat),]

# There is no line-key present. Thus we must do a manual line-calculation, based on stop dates.
other_CT = tolower(Hmisc::Cs(Cisplatin, Cyklofosfamid, Doxorubicin, Estramustin, 
                             Etoposid, Gemcitabin, Karboplatin, Metotrexat, Mitoxantron))
m1crpc_drugs <- c(Hmisc::Cs(abirateron, enzalutamid, "radium(Ra-223)diklorid", kabazitaxel, docetaxel), other_CT)
m1crpc_df <- drug_df %>% filter(drug_name %in% m1crpc_drugs)
other_drugs_df <- drug_df %>% filter(!drug_name %in% m1crpc_drugs)
m1crpc_df$other_CT = NA
m1crpc_df$other_CT[m1crpc_df$drug_name %in% other_CT] = T

dst_base <- haven::read_sas(paste0("R:\\PCBase41\\Works41\\PPC\\2018\\_Xofigostudie\\Data_source\\PPC\\sos_ppc1_cleaned_dst_ba_lopnr.sas7bdat")) # base

m1crpc_df[m1crpc_df==""]=NA


# We want Uts?ttning to be the last row in each line. event_dat could be the same for ins?ttning and uts?ttning, 
# and as ? is after U, we rename '?terins?ttning' before using arrange.
m1crpc_df$event <- gsub("?terins?ttning", "Aterins?ttning", m1crpc_df$event)

### 
# Throw away uncleanable data
# Identify men with start dates after stop dates, or remove one of duplicate stop dates (last one)

merge(m1crpc_df, dst_met[,c("LopNr","firstmet_all")], all.x=T) -> m1crpc_df 

### Calculate line ----------------------------------------
# Outer ldply splits on individual, inner splits on drug_name. 
m1crpc_df$death_dat <- as.Date(m1crpc_df$death_dat)
m1crpc_df$event_dat <- as.Date(m1crpc_df$event_dat)
m1crpc_df_s <- split(m1crpc_df, factor(m1crpc_df$LopNr))

# lpnr 556333, 427294
# for(i in 1:length(m1crpc_df_s)){if(any(m1crpc_df_s[[i]]$LopNr %in% "788801") ) break}
# 1657
bogus <- c(459055,165746,580226,16842,264227,88704,911162,556333,427294)

ldply(m1crpc_df_s, .progress = "text", .id=NULL, function(ind){
  tryCatch({
    # ind <- m1crpc_df_s[[1657]]
    ind$erase <- FALSE
    
    xx <- split(ind, with(ind, factor(paste(LopNr, drug_name))))
    
    ldply(xx, function(x){
      # x <- split(ind, factor(ind$drug_name!="docetaxel"))[[1]]
      
      #########
      # Correct a few incorrectly entered Xofigo-treatments have uts?ttning instead of paus. 
      # We set all but the last occurring uts?ttning to pause, if there are at least two 
      # uts?ttning and all of the stop dates lie within 60 days. 
      #########
      if(x$drug_name[1] %in% "radium(Ra-223)diklorid" & sum(x$event == "Uts?ttning")>2) {
        
        x %<>% arrange(event_dat)
        if(all(diff.Date(x$event_dat[x$event=="Uts?ttning"])<60)){
          x$event[x$event=="Uts?ttning"][-which.max(x$event_dat[x$event=="Uts?ttning"])] = "Paus"
        }
      }
      
      ######### 
      ### Swap index-explanation
      # Ordering events of the same date
      # Note the following two scenarios:
      # A) ("One day line of treatment", e.g. docetaxel one cycle)
      # 2000-01-01- Ins?ttning
      # 2000-01-01- Uts?ttning
      
      # B) ("Stop, but start again the same day", common for incorrectly entered Xofigo, and dose adjustment, etc)
      # 2000-01-01- Uts?ttning
      # 2000-01-01- Ins?ttning
      
      # To distinguish these from another (i.e. in what order the two events should be) we search backwards
      # for the closest event. If this is Uts?ttning, we treat the event as an A), otherwise as scenario B).
      # Technically - if there are two "uts?ttning" following eachother, we try to insert any non-uts?ttning-event
      # of the same date between the two "us?ttning". This is what "swap_index" does.
      #########
      
      x %<>% arrange(event_dat) %>% mutate(swap_index=NA)
      swap_these <- str_locate(paste0(as.numeric(x$event == "Uts?ttning"), collapse=""), "11")[1]
      x$swap_index[swap_these] = swap_these
      x %<>% arrange(event_dat, swap_index, event)
      
      x$line_count = 0
      # The first row is always the start of the first line
      x$line_count[1] = 1
      # Next, lines with "Uts?ttning" will have a new line on the NEXT succeeding row.
      # Note that the last row may or may not be Uts?ttning: If it is, any succeeding row would be a new line,
      # but there is no succeeding row. If the last line is not "Uts?ttning", there is no new line. 
      # I.e. the last line never matters, hence [-nrow(x)] below. 
      
      # If there is "Uts?ttning" without prior "Ins?ttning" then erase the whole person (i.e. all drug data)
      # Only look at if first uts?ttning is also first row in x
      if(sum(x$event %in% "Uts?ttning")>0){
        indx <- min(which(x$event %in% "Uts?ttning"))
        if(indx==1){x$erase <- TRUE}
        
      }
      
      x$line_count[which(x$event[-nrow(x)] %in% "Uts?ttning") + 1] = 1
      
      return(x)
    }, .id=NULL) -> ind_df
    
    
    
    ind_df$split_ind = NA
    ind_df$split_ind <- cumsum(ind_df$line_count)
    
    # Now the grouping of lines are found, but not the order of the lines. For this, we define "line_start" as the first date in each line.
    ind_df %<>% group_by(split_ind) %>% mutate(line_start = min(event_dat)) %>% ungroup()
    
    # Now, order the lines according to the line_start and calculate the line with cumsum. 
    ind_df %<>% arrange(line_start, drug_name, event_dat, swap_index, event) %>% select(LopNr, everything(), -swap_index) %>% as.data.frame()
    ind_df %<>% mutate(line = cumsum(line_count))
    
   
    
    ###### other_CT with preceding line also other_CT is not a new line
    # We build a variable other_CT_false_line which will keep how many extra lines that the other CT has created,
    # and in the end we remove this variable from line to get the correct line.
    # The initial NA pushes the 1:s in place of first line of new lines
    ind_df$other_CT_false_line[ind_df$other_CT %in% TRUE] <- c(NA, diff(ind_df$line[ind_df$other_CT %in% TRUE]))
    ind_df$other_CT_false_line[is.na(ind_df$other_CT_false_line)]=0
    ind_df$other_CT_false_line=cumsum(ind_df$other_CT_false_line) 
    ind_df$line <- ind_df$line - ind_df$other_CT_false_line
    
    # Cant have different line start dates for the same line:
    if(any(!is.na(ind_df$other_CT)>0 & !is.na(ind_df$other_CT_false_line))){
      ind_df[!is.na(ind_df$other_CT),] %<>% group_by(line) %>% mutate(line_start=min(c(event_dat,Inf), na.rm=T)) %>% 
        ungroup()
    }
    ###### 
    
   
    # Clean up
    ind_df %<>% select(-split_ind, -line_count)
    ind_df$error_message=NA
    
    ### Remove first line docetaxel if initiated within nine month of pre-gnrh-metastasis (NPCR or PPC)
    if(nrow(ind_df %>% filter(drug_name %in% "docetaxel" & line %in% 1) %>%
            filter(abs(line_start - firstmet_all) <= 9*30.5 & firstmet_all < line_start &
                   !is.na(firstmet_all) & !is.na(line_start))) > 0){
      
      # Set such docetaxel to line NA
      ind_df[(ind_df$drug_name %in% "docetaxel" & ind_df$lin %in% 1),"line"] = NA
      
      # 'These' then need line adjustment by -1
      these <- !(ind_df$drug_name %in% "docetaxel" & is.na(ind_df$line))
      ind_df$line[these] <- ind_df$line[these] - 1
    }
    
  
    # # Add line_cens=T/F and line_stop=date (line_cens=F for those with uts?ttning or death, 
    # # while line_cens=T means not dead and no uts?ttning)
     ind_df$line_cens = FALSE
     ind_df$line_stop = NA
     
     x_list <- split(ind_df, factor(paste0(ind_df$LopNr, ind_df$drug_name, ind_df$line)))
     
     ind_df <- ldply(x_list, function(x){
       # x <- x_list
       
       if( sum(x$event %in% "Uts?ttning" | !is.na(x$death_dat) )>0 ){
         yy <- c(x$event_dat[x$event %in% "Uts?ttning"], x$death_dat)
         if(sum(is.na(yy))==length(yy)){
           x$line_stop = as.Date("2018-12-31")
           x$line_cens = TRUE 
         } else{
           x$line_stop <- yy[which.min(as.numeric(yy))][1]
         }
      
       } else {
         x$line_stop = as.Date("2018-12-31")
         x$line_cens = TRUE 
       }
       return(x) 
     }) 
    


    if(sum(ind_df$erase)>0){ind_df$erase <- TRUE}
    ind_df
  }, error=function(e, ind_df){data.frame(error_message=as.character(e))})
}) -> line_df

ers <- unique(line_df$LopNr[line_df$erase])

line_df <- line_df[!line_df$erase,]

# line_df2 <- line_df[! line_df$drug_name %in% "Other CT",]
# sort(table(line_df2$LopNr),decreasing = TRUE)[1:10]


line_df %<>% arrange(LopNr, line_start, line)


if(all(is.na(line_df$error_message))){
  line_df %<>% select(-error_message)
} else {
  LopNr_with_errors <- line_df$LopNr[which(!is.na(line_df$error_message))]
  stop("Errors found in line calculation")
}

### Quality Control ----------------------------------------  
# As there was no line key, we make all quality checks I could think of. Is no better solution available.

# The removal of docetaxel first line was succesful
all(is.na(line_df[line_df$drug_name=="docetaxel" & line_df$line %in% NA & abs(ymd(line_df$line_start) - ymd(line_df$firstmet_all)) < 9*30.5 & 
                    line_df$firstmet_all < line_df$line_start,"line"]))

line_df$event_dat <- as.character(line_df$event_dat)
line_df$death_dat <- as.character(line_df$death_dat)


line_df <- bind_rows(line_df, other_drugs_df)


# Quality control - how many of each line did we get, compared to the previous. 
dst_drug <-  haven::read_sas(paste0("R:\\PCBase41\\Works41\\PPC\\2018\\_Xofigostudie\\Data_source\\PPC\\sos_ppc4_cleaned_dst_dr_lopnr.sas7bdat")) # load("F:\\PPC\\PPC_41\\Convenient\\rdata_versioner\\dst_drug.rdata")
line_df[line_df$drug_name %in% other_CT, "drug_name"] = "Other CT"


new <- line_df %>% filter(LopNr %in% dst_drug$LopNr) %>% distinct(LopNr, drug_name, line) %>% select(-LopNr) %>%  table()

dst_drug$drug_name <- tolower(dst_drug$drug_name)

length(unique(line_df$LopNr))
length(unique(dst_drug$LopNr))

dst_drug[dst_drug$drug_name %in% other_CT, "drug_name"] = "Other CT"
# dst_drug-file was extracted from PPC at a later date due to a data error. We remove the late entries
# dst_drug <- dst_drug %>% filter(drug_name %in% c("other ct", "ra223", unique(line_df$drug_name))) %>% 
#   filter(drugstart_dat <= as.Date("2018-01-26"))

old <- dst_drug %>% filter(drug_name %in% c("other ct", "ra223", unique(line_df$drug_name))) %>%  distinct(LopNr, drug_name, drug_line_crpc) %>% select(-LopNr) %>%  table()
old <- as.data.frame.matrix(old[rowSums(old)!=0,])
new <- as.data.frame.matrix(new[rowSums(new)!=0,])
old <- old[rownames(old) %in% c(rownames(new),"ra223","other ct"),]
# Close enough. 
cbind(old, rowSums(old))
cbind(new, rowSums(new))
sum(old)
sum(new)

drug_df <- line_df %>% select(-.id, -firstmet_all, -other_CT_false_line, -other_CT)


cleaded_xof <- unique(dst_drug$LopNr[dst_drug$drug_name %in%  "ra223"])
length(cleaded_xof)

new_xof <- unique(drug_df$LopNr[drug_df$drug_name %in% "radium(Ra-223)diklorid"])
length(new_xof)

sum(!new_xof %in% cleaded_xof)
odd_xof <- new_xof[!new_xof %in% cleaded_xof]
# gnrh_df[gnrh_df$LopNr %in% odd_xof,]


drug_df[drug_df$LopNr %in% odd_xof[1],]
dst_drug[dst_drug$LopNr %in%  odd_xof[1],]
m1crpc_df[m1crpc_df$LopNr %in% odd_xof[1], ]
drugs[drugs$LopNr %in% odd_xof[1],c("sub_substans_vd_sub_subnamnrek_","sub_substans_vd_sub_subnamnrek","LMIH_handelseorsak_Beskrivning","LMIH_handelse_Beskrivning","LMIH_handelsedatum","LMIH_handelseorsakann",
                                    "LMIU_utsattningsorsak_Beskrivni","LMIU_utsattningsdatum","LMIU_utsattningsorsakann")]

drug_df[drug_df$LopNr %in% odd_xof[2],]
dst_drug[dst_drug$LopNr %in%  odd_xof[2],]
m1crpc_df[m1crpc_df$LopNr %in% odd_xof[2], ]
drugs[drugs$LopNr %in% odd_xof[2],c("sub_substans_vd_sub_subnamnrek_","sub_substans_vd_sub_subnamnrek","LMIH_handelseorsak_Beskrivning","LMIH_handelse_Beskrivning","LMIH_handelsedatum","LMIH_handelseorsakann",
                                    "LMIU_utsattningsorsak_Beskrivni","LMIU_utsattningsdatum","LMIU_utsattningsorsakann")]

## Save ----------------

save(drug_df, file="R:\\PCBase41\\Works41\\PPC\\2018\\_Xofigostudie\\Data_source\\PPC_fixed_drug\\cleaned_drugs_df.Rdata")







#