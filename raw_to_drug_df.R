###
# By Oskar Gauffin and modified by Marcus Westerberg
# Modified: 2019-08-09

# !diagnostics off
# The line above prevents dplyr from outputting a lot of nonsense warnings

##############################################
# Line calculation and drug duration-script
##############################################
# This script is written to calculate crpc-line of treatments, and from this, treatment duration, 
# based on the cessation (uts?ttning) of the drug or death date of the patient.

# The overall improvement compared to the previous dst_drug-file (line was from first date to last date of the drug) is to find multiple lines for the same drug
# this was not the case in the previous script. That is, one day docetaxel, ten year abiraterone, one day docetaxel should be considered as three lines,
# not as a first line of docetaxel with duration 10 years and 2 days, and an overlapping "second" line of abiraterone.

### A few special cases are built in:
# Docetaxel started within nine month of pre-GnRH-metastases (either from NPCR (at dx) or PPC (later)) in the first line, is not counted
# as a crpc-treatment.

# Other CT (CM, KEES, etc) could be multi-drug-line. Thus, secondary initiations of other CT only counts as a new line of treatment if some other 
# drug (not in "Other CT") is in between two initiations of "Other CT"-drugs. 

# Xofigo is sometimes mistakenly ended instead of paused after each cycle. We try to fix this. 

# Initiation and cessation could occur at the same date, and then preceding cessations will decide
# whether this is a one-day-line, or an ended line and a new line started on the same date. 

# It was not considered sensible to require that the initiation is marked as an initiation, as some have free text comments or "side effects" as 
# initiation comment (probably, the initiation was due to some side effect of another drug, or possible entered retrospectively that the drug was initiated
# but with side effects found later.)

### Convert raw-PPC to usable PPC-INCA-files.
### Data cleaning -----------------------------
packages <- c("magrittr","lubridate","plyr","dplyr", "stringr","data.table")
DIV::library()

# Ladda l?kemedelsvyn (allt i en df, dvs efter mergning som g?rs i vy-steget)
load("C:/Users/gao002/Desktop/data_till_drug_df/PPC_PCBase_extraction_drug.RData")
drugs <- df
load("C:/Users/gao002/Desktop/data_till_drug_df/zzz_R_root_persnr_kopia.RData")
dst_pers <- df
dst_pers %<>% rename(lpnr=U_pk) %>% select(lpnr, PERSNR)
drugs <- merge(drugs, dst_pers, by="PERSNR", all.x=T)

# Load dst_base and dst_diag
dat <- gsub("-","", Sys.Date())
load(paste0("Z:/Uppsala/oc2foya/ppc/temp/PPC_dst_base_", dat, ".rdata"))
load(paste0("Z:/Uppsala/oc2foya/ppc/temp/PPC_dst_diag_", dat, ".rdata"))

dst_base %<>% select(lpnr, birth_dat, death_dat, lpnr, firstmet_dat) %>% rename(firstmet_dat_base = firstmet_dat)

###############################
# Metastases from PPC
###############################
# For PPC we only care about metastases occurring before the initiation of GnRH.
# The latest of those metastases could explain initiation of docetaxel, if it is within nine months before docetaxel initiation.
# Create df with metastases, dst_met, out of dst_base:
dst_met <- dst_base %>% select(lpnr, firstmet_dat_base) %>% filter(!is.na(firstmet_dat_base)) %>% 
  select(lpnr, firstmet_dat_base) %>% rename(met_dat = firstmet_dat_base)

# Set missing data to NA
dst_diag$image_metast[dst_diag$image_metast %in% ""] <- NA
dst_diag$created_dat[dst_diag$created_dat %in% ""] <- NA

# Select the metastases in dst_diag.
dst_diag %<>% filter(image_assess %in% c("Blandad Respons", "Progress", "1:a metastas") | !is.na(image_metast)) %>% 
  select(lpnr, image_dat) %>% rename(met_dat=image_dat)

######################
# Find first initiation of GnRH in drugs
######################
gnrh <- tolower(c("Leuprorelin","Goserelin","Buserelin","Triptorelin","Histrelin"))
drugs$LMIH_handelsedatum[drugs$LMIH_handelsedatum %in% "NA"] <- NA
drugs$LMIH_handelsedatum <- as.Date(drugs$LMIH_handelsedatum)

# Marcus added "Inf" in the min, to avoid a few warnings when people had no first gnrh.
gnrh_df <- drugs %>% rename(drug_name = sub_substans_vd_sub_subnamnrek_1) %>% 
  filter(drug_name %in% gnrh) %>% group_by(lpnr) %>% 
  mutate(first_gnrh = min(c(LMIH_handelsedatum, Inf), na.rm=T)) %>% 
  select(lpnr, first_gnrh)

# Those with Inf are then removed in the next step.
if(sum(gnrh_df$first_gnrh==Inf)!=0){
  gnrh_df <- gnrh_df[!gnrh_df$first_gnrh==Inf,]
}


# Merge dst_met with metastases from dst_diag, and first gnrh from gnrh_df:
dst_met <- rbind.data.frame(dst_met, dst_diag) %>% distinct()
dst_met <- merge(dst_met, gnrh_df, all.x=T)

# Mark metastases found after initiation of GnRH:
met_after_gnrh <- !is.na(dst_met$first_gnrh) & !is.na(dst_met$met_dat) &
  dst_met$first_gnrh <= dst_met$met_dat

dst_met$post_GnRH = FALSE
dst_met$post_GnRH[met_after_gnrh] = TRUE

###############################
# Metastases from NPCR
###############################
# For NPCR we use all metastases (as no one starts gnrh before date of diagnosis)
load("C:/Users/gao002/Desktop/data_till_drug_df/PPC_M1_info_from_NPCR.RData")
npcr <- df
npcr <- merge(npcr, dst_pers, by="PERSNR")
npcr$d0_m_txt <- factor(npcr$d0_m, labels=c("M0","M1","M0"))
npcr %<>% mutate_if(is.factor, as.character) 
npcr %<>% select(lpnr, d0_date, d0_m_txt) %>% filter(d0_m_txt=="M1") %>% select(lpnr, met_dat=d0_date, -d0_m_txt) %>% 
  mutate(met_dat = substr(met_dat, 0, 10))
npcr$post_GnRH = F
npcr$first_gnrh=NA

dst_met <- rbind.data.frame(dst_met, npcr)
dst_met$met_dat <- as.Date(dst_met$met_dat)

dst_met %<>% distinct()
dst_met$met_dat[dst_met$met_dat %in% "NA"] = NA
dst_met <- dst_met[!is.na(dst_met$met_dat),]

load("C:/Users/gao002/Desktop/data_till_drug_df/PPC_PCBase_extraction_drug.RData")
drugs <- df
drugs <- merge(drugs, dst_pers, by="PERSNR", all.x=T)
drugs <- merge(dst_base, drugs, all.y=T)
drugs <- drugs %>% rename(drug_name = sub_substans_vd_sub_subnamnrek_1)
drugs$drug_name <- gsub(" \\(vattenfri\\)", "", drugs$drug_name)

start_var <- c("lpnr","drug_name", "LMIH_handelseorsak_Beskrivning","LMIH_handelse_Beskrivning",
               "LMIH_handelsedatum", "LMIH_handelseorsakann", "death_dat")
stop_var <- c("lpnr","drug_name", "LMIU_utsattningsdatum", "LMIU_utsattningsorsak_Beskrivning",
              "LMIU_utsattningsorsakann", "death_dat")

drug_df <- bind_rows(drugs[, start_var] %>% rename(event = LMIH_handelseorsak_Beskrivning,
                                                   comment=LMIH_handelse_Beskrivning,
                                                   event_dat=LMIH_handelsedatum,
                                                   free_text_comment=LMIH_handelseorsakann) %>% 
                       mutate_if(is.factor, as.character),
                     
                     drugs[, stop_var] %>% rename(comment=LMIU_utsattningsorsak_Beskrivning,
                                                  event_dat=LMIU_utsattningsdatum,
                                                  free_text_comment=LMIU_utsattningsorsakann)) %>%
  filter(!is.na(event_dat)) %>% 
  mutate_if(is.factor, as.character)

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
m1crpc_df[m1crpc_df==""]=NA

# Load dst_base
load(paste0("Z:/Uppsala/oc2foya/ppc/temp/PPC_dst_base_", dat,".rdata"))

# We want Uts?ttning to be the last row in each line. event_dat could be the same for ins?ttning and uts?ttning, 
# and as ? is after U, we rename '?terins?ttning' before using arrange.
m1crpc_df$event <- gsub("?terins?ttning", "Aterins?ttning", m1crpc_df$event)

### Calculate line ----------------------------------------
# Outer ldply splits on individual, inner splits on drug_name. 
m1crpc_df$death_dat <- as.Date(m1crpc_df$death_dat)
m1crpc_df$event_dat <- as.Date(m1crpc_df$event_dat)
m1crpc_df_s <- split(m1crpc_df, factor(m1crpc_df$lpnr))

line_df <- ldply(m1crpc_df_s, .progress = "text", .id=NULL, function(ind){
  
  # tryCatch stores the error messages, and makes ldply run through all records.
  tryCatch({
    # ind <- m1crpc_df_s[[i]]
    ind$erase <- FALSE
    ind$error_message <- NA
    
    drug_ind <- split(ind, with(ind, factor(paste(lpnr, drug_name))))
    
    ldply(drug_ind, function(x){
      # x <- drug_ind[[1]]
      
      #########
      # Correct a few incorrectly entered Xofigo-treatments with uts?ttning instead of paus. 
      # We set all but the last occurring uts?ttning to pause, if there are at least two 
      # uts?ttning and all of the stop dates lie within 60 days.
      # This rule of thumb is from Ingela, as it's not likely that patients have more than two 
      # lines of Xofigo.
      #########
      if(x$drug_name[1] %in% "radium(Ra-223)diklorid" & sum(x$event %in% c("Uts?ttning","Utsättning"))>2) {
        
        x %<>% arrange(event_dat)
        if(all(diff.Date(x$event_dat[x$event %in% c("Uts?ttning","Utsättning")])<60)){
          x$event[x$event %in% c("Uts?ttning","Utsättning")][-which.max(x$event_dat[x$event %in% c("Uts?ttning","Utsättning")])] = "Paus"
        }
      }
      
      ######### 
      ### Swap index-explanation
      # Ordering events of the same date
      # Note the following two scenarios:
      # A) ("One day line of treatment", e.g. docetaxel one cycle with cessation due to side effects)
      # 2000-01-01- Ins?ttning
      # 2000-01-01- Uts?ttning
      
      # B) ("Stop, but start again the same day", common for incorrectly entered Xofigo, etc)
      # 2000-01-01- Uts?ttning
      # 2000-01-01- Ins?ttning
      
      # To distinguish these from another (i.e. in what order the two events should be) we search backwards
      # for the closest event. If this is Uts?ttning, we treat the event as an A), otherwise as scenario B).
      # Technically - if there are two "uts?ttning" following eachother, we try to insert any non-uts?ttning-event
      # of the same date between the two "us?ttning". This is what "swap_index" does.
      #########
      
      x %<>% arrange(event_dat) %>% mutate(swap_index=NA)
      swap_these <- str_locate(paste0(as.numeric(x$event %in% c("Uts?ttning","Utsättning")), collapse=""), "11")[1]
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
      if(sum(x$event %in% c("Uts?ttning","Utsättning"))>0){
        indx <- min(which(x$event %in% c("Uts?ttning","Utsättning")))
        if(indx==1){x$erase <- TRUE}
        
      }
      
      x$line_count[which(x$event[-nrow(x)] %in% c("Uts?ttning","Utsättning")) + 1] = 1
      
      return(x)
    }, .id=NULL) -> ind_df
    
    ind_df$split_ind = NA
    ind_df$split_ind <- cumsum(ind_df$line_count)
    
    
    # find dates that are not "utsättning"
    ind_df$index <- 1:nrow(ind_df)
    ind_df2 <- ind_df[ ! ind_df$event %in% c("Uts?ttning","Utsättning"), ]
    
    # Now the grouping of lines are found, but not the order of the lines. For this, we define "line_start" as the first date in each line, excluding "utsättning".
    ind_df2 %<>% group_by(split_ind) %>% mutate(line_start = min(event_dat))  %>% ungroup()
    ind_df2 <- ind_df2[,c("index","line_start")]
    ind_df <- merge(ind_df,ind_df2,by="index",all.x=TRUE)
    
    ind_df %<>% group_by(split_ind) %>% mutate(line_start = min(line_start,na.rm = TRUE))  %>% ungroup()
    
    # erase if line_start is missing
    if(sum(is.na(ind_df$line_start))>0 ){ ind_df$erase = TRUE }
    
    # Now, order the lines according to the line_start and calculate the line with cumsum. 
    ind_df %<>% arrange(line_start, drug_name, event_dat, swap_index, event) %>% select(lpnr, everything(), -swap_index) %>% as.data.frame()
    ind_df %<>% mutate(line = cumsum(line_count))
    
    ###### other_CT with preceding line that is also other_CT is not a new line.
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
    
    # Clean up
    ind_df %<>% select(-split_ind, -line_count)
    # ind_df$error_message=NA
    
    ### Adjust first line docetaxel to NA, if initiated within nine month of pre-gnrh-metastasis (NPCR or PPC) 
    ### and abscence of metastases from GnRH-initiation to docetaxel-start, and docetaxel-start after gnrh-start.
    if(any(ind_df$drug_name %in% "docetaxel" & ind_df$line %in% 1)){
      
      ind_doc_line1 <- ind_df$line_start[ind_df$drug_name %in% "docetaxel" & ind_df$line %in% 1][1]
      dst_met_ind <- dst_met[dst_met$lpnr == ind_df$lpnr[1],]
      
      
      if(
        # if there are any pre-gnrh-metastases
        nrow(dst_met_ind[!dst_met_ind$post_GnRH,]) > 0 & 
        
        # and if there are pre-gnrh metastases within 9 months of docetaxel initiation
        any(abs(dst_met_ind$met_dat[!dst_met_ind$post_GnRH] - ind_doc_line1) < 9*30.5) &
        
        # and not any post-gnrh, pre-docetaxel metastases,
        ! any(dst_met_ind$met_dat[dst_met_ind$post_GnRH] < ind_doc_line1) &
        
        # and docetaxel is initiated after gnrh, or there is no gnrh at all.
        (dst_met_ind$first_gnrh[1]  < ind_doc_line1 | all(is.na(dst_met_ind$first_gnrh))) |
        
        # or if docetaxel is initiated before gnrh - then it's certainly not crpc.
        (ind_doc_line1 <= min(c(dst_met_ind$first_gnrh[1])) )
      ){
        
        # Set such docetaxel to line NA
        ind_df[(ind_df$drug_name %in% "docetaxel" & ind_df$lin %in% 1),"line"] = NA
        
        # Succeeding lines ('These' below) then need line adjustment by -1
        these <- !(ind_df$drug_name %in% "docetaxel" & is.na(ind_df$line))
        ind_df$line[these] <- ind_df$line[these] - 1
      }
    }
    
    # # Add line_cens=T/F and line_stop=Sys.Date (line_cens=F for those with uts?ttning or death, 
    # # while line_cens=T means not dead and no uts?ttning)
    ind_df$line_cens = FALSE
    ind_df$line_stop = Sys.Date()
    
    x_list <- split(ind_df, factor(paste0(ind_df$lpnr, ind_df$drug_name, ind_df$line)))
    
    ind_df <- ldply(x_list, function(x){
      # x <- x_list
      
      # If there is a cessation or death:
      if(sum(x$event %in% c("Uts?ttning","Utsättning") | !is.na(x$death_dat)) > 0){
        
        yy <- c(x$event_dat[x$event %in% c("Uts?ttning","Utsättning")], x$death_dat)
        
        # Cessation might be without date. Check if cessation date and death dates are all NA. 
        # In that case, set censoring = T at current date.
        if(sum(is.na(yy))==length(yy)){
          x$line_stop = Sys.Date()
          x$line_cens = TRUE 
        } else {
          # If there is cessation or death, we store the earliest one and do not change line_cens from FALSE:
          x$line_stop <- yy[which.min(as.numeric(yy))][1]
        }
        
        # If there were no cessation and no deaths, set censoring = T and use todays date:
      } else {
        x$line_stop = Sys.Date()
        x$line_cens = TRUE 
      }
      return(x) 
    })
    
    if(sum(ind_df$erase)>0){
      ind_df$erase <- TRUE
    }
    
    ind_df
  }, error=function(e, ind_df){data.frame(error_message=as.character(e))})
})

# unique(m1crpc_df$lpnr[!m1crpc_df$lpnr %in% line_df$lpnr])
table(line_df$erase,useNA = "always")

# Patient with cessation without preceding initiation are excluded - not meaningful to try to fix those. 
ers <- unique(line_df$lpnr[line_df$erase])
if(sum(is.na(ers))>0){print("Some missing lpnr!");print(sum(is.na(ers)))}

line_df <- line_df[!line_df$erase,]

# Remove NA-rows (ok since these have not correct data)
line_df <- line_df[!is.na(line_df$lpnr),]

# Patient with cessation without preceding initiation are excluded - not meaningful to try to fix those. 
ers <- unique(line_df$lpnr[line_df$erase])
line_df <- line_df[!line_df$erase,]
line_df %<>% arrange(lpnr, line, event_dat)

# This is for debugging. While running on INCA I would prefer a running script before a script that breaks down at the first error.
# if(all(is.na(line_df$error_message))){
#   line_df %<>% select(-error_message)
# } else {
#   lpnr_with_errors <- line_df$lpnr[which(!is.na(line_df$error_message))]
#   stop("Errors found in line calculation")
# }

line_df$event_dat <- as.character(line_df$event_dat)
line_df$death_dat <- as.character(line_df$death_dat)

# Put back all non-CRPC drugs (bika, GnRH and so on)
line_df <- bind_rows(line_df, other_drugs_df)
line_df[line_df$drug_name %in% other_CT, "drug_name"] = "Other CT"

# Remove demo-patients and NA-rows
line_df <- merge(line_df, dst_pers, by="lpnr", all.x=T)
line_df <- line_df[substr(line_df$PERSNR, 0, 2)  %in% c("19","20"),]
line_df <- line_df[!is.na(line_df$lpnr),]

## Save ----------------

colnames(line_df)
drug_df <- line_df %>% select(-.id, -erase, -other_CT_false_line, -other_CT, -PERSNR)
save(drug_df, file="D:/Rskript/Uppsala/oc2foya/PPC/temp/drug_df.Rdata")

############################################# 
# Compare line_df (new) to dst_drug (old)
#############################################
# new <- line_df %>% filter(lpnr %in% dst_drug$lpnr) %>% distinct(lpnr, drug_name, line) %>% select(-lpnr) %>%  table()
# new <- as.data.frame.matrix(new[rowSums(new)!=0,])
# 
# dat <- gsub("-","",Sys.Date())
# load(paste0("Z:/Uppsala/oc2foya/ppc/temp/PPC_dst_drug_",dat,".rdata"))
# dst_drug$drug_name <- tolower(dst_drug$drug_name)
# docetaxel_line_0 <- dst_drug %>% filter(drug_name == "docetaxel" & is.na(drug_line_crpc)) %>% distinct(lpnr) %>% nrow()
# dst_drug[dst_drug$drug_name %in% other_CT, "drug_name"] = "Other CT"
# old <- dst_drug %>% filter(drug_name %in% c("other ct", "ra223", unique(line_df$drug_name))) %>%  distinct(lpnr, drug_name, drug_line_crpc) %>% select(-lpnr) %>%  table()
# old <- as.data.frame.matrix(old[rowSums(old)!=0,])
# old <- old[rownames(old) %in% c(rownames(new),"ra223","other ct"),]
# old <- cbind("0"=c(0,docetaxel_line_0, rep(0,4)), old)
# 
# new
# old
# new[,2:8] - old[,2:8]
# # Sv?rt att definiera vad som ?r "en" ?ndring. Om f?rsta linjen ?ndras f?r n?gon med tre linjer => ?ndring p? fyra platser i tabellen.
# 
# 
# # Confusing that more patients have docetaxel in first line in new script. Check manually:
# new_abi_1 <- line_df %>% filter(line==1 & drug_name %in% "docetaxel") %>% select(lpnr) %>% distinct(lpnr)
# old_abi_1 <- dst_drug %>% filter(drug_line_crpc==1 & drug_name == "docetaxel") %>% select(lpnr) %>% distinct()
# check_these <- new_abi_1[!new_abi_1$lpnr %in% old_abi_1$lpnr,]
# 
# check_df <- data.frame("lpnr"=check_these)
# check_df$cause = NA
# 
# # ind_lpnr = XXXX
# dst_drug[dst_drug$lpnr %in% ind_lpnr & dst_drug$drug_name %in% "docetaxel","drug_line_crpc"]
# line_df[line_df$lpnr %in% ind_lpnr & line_df$drug_name %in% "docetaxel","line"]
# 
# for(i in 1:nrow(check_df)){
#   # i = 98
# 
#   # 221 out of 297 had intermittent metastases. 
#   ind_lpnr=check_df$lpnr[i]
#   (gnrh_dat <- min(gnrh_df[gnrh_df$lpnr==ind_lpnr,]$first_gnrh, na.rm=T))
#   (dst_met_ind <- dst_met[dst_met$lpnr==ind_lpnr,])
#   (doc_dat <- line_df[line_df$lpnr %in% ind_lpnr & line_df$drug_name %in% "docetaxel", "line_start"])
#   
#   if(any((min(gnrh_dat) < dst_met_ind$met_dat) & (dst_met_ind$met_dat < min(doc_dat)))){
#     check_df$cause[i] <- "int_met"
#   }
#   
#   # Check if line has been splitted compared to dst_drug. 64 out of those not above (297-221) had splitted lines.
#   if(is.na(check_df$cause[i]) & length(dst_drug[dst_drug$lpnr %in% ind_lpnr & dst_drug$drug_name %in% "docetaxel","drug_line_crpc"]) < 
#      length(line_df[line_df$lpnr %in% ind_lpnr & line_df$drug_name %in% "docetaxel","line"])){
#     check_df$cause[i] = "splitted line"
#   }
# }
# 
# # 7 lpnr do not fall into these two categories. We leave them as they are. 
# check_df[is.na(check_df$cause),"lpnr"]
# table(check_df$cause)