# November 2022
# Mónica Leiria de Mendonça
# https://github.com/mlmms

# search for "OUTPUT LIST" to go to relevant outputs for this analysis
library(survival)
library(survminer)
library(stringr)
library(glmSparseNet)
set.seed = 1000

#SINGLE OMICS
#SINGLE OMICS
#SINGLE OMICS
#SINGLE OMICS
#SINGLE OMICS

# DATA LOADING -----------------------------------------------------------------

#////////////////////////////////////////////////////////////////////////
#///////////////////// LGG + GBM / PANGLIOMA analysis////////////////////
#////////////////////////////////////////////////////////////////////////

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~ lgg+gbm ALL genes analysis ~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2016
load("~/paraRoberta/single-omics-ds/single-omics_final_ds_2016.RData")            #643 p X 802 var
Panglioma_2016_ALL <- Pan_glioma_f
rm(Pan_glioma_f, Pan_glioma_norm_f)

#2021
load("~/paraRoberta/single-omics-ds/recent_DS/single-omics_final_ds_2021.RData")  #619 p X 1002 var
Panglioma_2021_ALL <- Pan_glioma_f
rm(Pan_glioma_f, Pan_glioma_norm_f)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~ lgg+gbm HUB genes analysis ~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2016
load("~/paraRoberta/single-omics-ds/recent_DS/single-omics_ex_sel_genes_2016.RData")  #643 p X 382 var
Panglioma_2016_HUB <- Pan_glioma_f
rm(Pan_glioma_f, Pan_glioma_norm_f)


#2021
load("~/paraRoberta/single-omics-ds/recent_DS/single-omics_ex_sel_genes_2021.RData")  #619 p X 551 var
Panglioma_2021_HUB <- Pan_glioma_f
rm(Pan_glioma_f, Pan_glioma_norm_f)


#///////////////////////////////////////////////////////
#///////////////////// LGG analysis ////////////////////
#///////////////////////////////////////////////////////

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~ lgg ALL genes analysis ~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2016
load("~/paraRoberta/single-omics-ds/recent_DS/single-omics-LGG_all_glasso2016.RData") # 494 p X 685 var
Lgg_2016_ALL <- LGG_all
rm(LGG_all, LGG_all_norm)

#2021
load("~/paraRoberta/single-omics-ds/recent_DS/single-omics-LGG_all_glasso2021.RData") # 420 p X 900 var
Lgg_2021_ALL <- LGG_all
rm(LGG_all, LGG_all_norm)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~ lgg HUB genes analysis ~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2016
load("~/paraRoberta/single-omics-ds/recent_DS/single-omics-LGG_ex_sel_glasso2016.RData") # 494 p X 325 var
Lgg_2016_HUB <- LGG_ex_sel
rm(LGG_ex_sel, LGG_ex_sel_norm)

#2021
load("~/paraRoberta/single-omics-ds/recent_DS/single-omics-LGG_ex_sel_glasso2021.RData") # 420 p X 531 var
Lgg_2021_HUB <- LGG_ex_sel
rm(LGG_ex_sel, LGG_ex_sel_norm)



#CLINICAL DATA
load("~/paraRoberta/clinical-Panglioma.RData") #panglioma_clinic_single_omics_2016 643 p X 16 var #panglioma_clinic_single_omics_2021 619 p X 13 var #rm(panglioma_clinic_single_omics_2016, panglioma_clinic_single_omics_2021, panglioma_clinic_multi_omics)
SURV.raw <- panglioma_clinic_single_omics_2016
rm(panglioma_clinic_single_omics_2016, panglioma_clinic_single_omics_2021, panglioma_clinic_multi_omics)


#SURVIVAL AND HISTOLOGICAL CLASSES DATA
#load("~/paraRoberta/clinical-dataset-glioma-subtypes_2021.RData")
#astro_clinic, oligo_clinic, gbm_clinic, unclassified_samples
#SURV.raw                    <- rbind(astro_clinic, oligo_clinic, gbm_clinic, unclassified_samples)
rownames(SURV.raw)          <- SURV.raw$Patient_ID
CLIN.raw                    <- SURV.raw
SURV.raw                    <- SURV.raw[,c(1,3,4,5)]               
#SURV_lgg.raw$patientID     <- rownames(SURV_lgg.raw) %>% toupper()
#SURV_lgg.raw$patientID    <- chartr('.', '-', SURV_lgg.raw$patientID)
#rownames(SURV_lgg.raw)    <- SURV_lgg.raw$patientID

SURV  <- dplyr::rowwise(SURV.raw) %>%
  dplyr::mutate(time = max(days_to_last_followup, days_to_death, na.rm = TRUE) ) %>%
  dplyr::select(Patient_ID, status = vital_status, time) # %>%    # Keep only survival variables and codes
#  dplyr::filter(!is.na(time) ) #removes 4 patients   # #& time > 10 Discard individuals with survival time less or equal to 10

SURV              <- as.data.frame(SURV)
rownames(SURV)    <- SURV$Patient_ID
SURV              <- SURV[,-1]
#rm(astro_clinic, gbm_clinic, oligo_clinic, unclassified_samples)

## MATCH and reorder SURV data (with more patients) according to the Panglioma data 2016 patients
Panglioma_idx   <- match(rownames(Panglioma_2016_ALL), rownames(SURV))
SURV_reordered  <- SURV[Panglioma_idx, ]

SURV_reordered$time <- as.numeric(SURV_reordered$time)
SURV_reordered$status <- as.numeric(SURV_reordered$status)

SURV <- SURV_reordered
rm(SURV_reordered)

SURV$Patient_ID <- rownames(SURV)
merged_SURV     <- merge(SURV, CLIN.raw[,c(1,11)], by.x = "Patient_ID")
#aux   <- match(rownames(SURV), rownames(SURV_with_histological_type))


# SURVIVAL ANALYSIS -----------------------------------------------------------
#----------------- PANGLIOMA 2016

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ KAPLAN MEYER MODEL ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FIT ONE KAPLAN MEIER CURVE FOR ALL CASES AS ONE GROUP "~1" #(based on my script survival_tcga_gbm.v1.R)
km <- with(SURV, Surv(time, status))
km_merged <- with(merged_SURV, Surv(time, status))

head(km,80) #"+" after the time indicates censoring.

# CENSORED data -------------
km.vector <- as.character(km)
length(grep("\\+", km.vector))

#count censored data by tumor type
km_merged_vector <- as.character(km_merged)
merged_SURV$km_merged_vector <- km_merged_vector

table(merged_SURV$histological_type)
gbm_SURV   <- merged_SURV[merged_SURV$histological_type == "untreated primary (de novo) gbm", ]
astro_SURV <- merged_SURV[merged_SURV$histological_type == "astrocytoma", ]
oligoastrocytoma_SURV <- merged_SURV[merged_SURV$histological_type == "oligoastrocytoma", ]
oligodendroglioma_SURV <- merged_SURV[merged_SURV$histological_type == "oligodendroglioma", ]

length(grep("\\+", gbm_SURV$km_merged_vector))
length(gbm_SURV$Patient_ID) 
pc_gbm <- 1/length(gbm_SURV$Patient_ID) * length(grep("\\+", gbm_SURV$km_merged_vector)) * 100 #% censored data
cat("% gbm censored cases = ", pc_gbm); rm(pc_gbm)

length(grep("\\+", astro_SURV$km_merged_vector))
length(astro_SURV$Patient_ID)
pc_astro <- 1/length(astro_SURV$Patient_ID) * length(grep("\\+", astro_SURV$km_merged_vector)) * 100 #% censored data
cat("% astrocytoma censored cases = ", pc_astro); rm(pc_astro)


length(grep("\\+", oligoastrocytoma_SURV$km_merged_vector))
length(oligoastrocytoma_SURV$Patient_ID)
pc_oligoastro <- 1/length(oligoastrocytoma_SURV$Patient_ID) * length(grep("\\+", oligoastrocytoma_SURV$km_merged_vector)) * 100 #% censored data
cat("% oligoastrocytoma censored cases = ", pc_oligoastro); rm(pc_oligoastro)


length(grep("\\+", oligodendroglioma_SURV$km_merged_vector))
length(oligodendroglioma_SURV$Patient_ID)
pc_oligodendro <- 1/length(oligodendroglioma_SURV$Patient_ID) * length(grep("\\+", oligodendroglioma_SURV$km_merged_vector)) * 100 #% censored data
cat("% oligodendroglioma censored cases = ", pc_oligodendro); rm(pc_oligodendro)



# (skip) FORMULA ------------------
#a formula object, which must have a Surv object as the response on the left
#of the ~ operator and, if desired, terms separated by + operators on the right.
#One of the terms may be a strata object. For a single survival curve the right
#hand side should be ~ 1.
#stan.KM <- survfit( Surv (time, status ) ~ 1 , data = SURV)
#plot(stan.KM, col=c(2,4,6), xlab = "Time (Days)", ylab = "Survival Probability (%)", main = 'Kaplan Meyer Survival Curve')  #col = colour
#autoplot(stan.KM, xlab = "Time (days)", ylab = "Survival Probability (%)", main = "Kaplan Meyer Survival curve of TCGA-PanGlioma cases_2016 (n=643)")

#summary <- summary(stan.KM, times = c(1,30,60,90*(1:10),1500,2000))
#summary
#capture.output(summary, file = "summary.txt")

# Change color, linetype by strata, risk.table color by strata
# ggsurvplot(stan.KM,
#            pval = TRUE, conf.int = TRUE,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            palette = c("#E7B800", "#2E9FDF"))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~ COX PROPORTIONAL HAZARDS MODEL ~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fit Cox Model

# #AUX covariates list - concatenate colnames vector into a string
# Pan_glioma_i_2016_2<-Pan_glioma_i_2016[complete.cases(Pan_glioma_i_2016),]
# #STRING_all_16k <- paste0(colnames(Pan_glioma_i_2016_2), collapse = " + ")
# 
# Panglioma_2016_ALL<-Pan_glioma_f_2016[complete.cases(Pan_glioma_f_2016),]
# #STRING_2016_643v <- paste0(colnames(Panglioma_2016_ALL), collapse = " + ")
# 
# #change "-" character in colnames to "aaa",
# #as it was giving trouble in the expression for Cox model fitting
# colnames(Panglioma_2016_ALL) <- gsub(x = colnames(Panglioma_2016_ALL), pattern = "\\-", replacement = "aaa")  
# #STRING_2016_643v <- paste0(colnames(Panglioma_2016_ALL), collapse = " + ")
# 
# Pan_glioma_i_2016_and_SURV <- cbind(Pan_glioma_i_2016_2, SURV)
# Pan_glioma_f_2016_and_SURV <- cbind(Panglioma_2016_ALL, SURV)
# 
# 
# # #NEW
# # #Compare:
# # # A.1 - Cox model with features selected from GLASSO (2016: 802 variables)
# # 
# # #g <- coxph(Surv(time, status) ~ A2BP1 + ABI3, data = Pan_glioma_f_2016_and_SURV)
# # glasso_cox_model_643v_2016 <- coxph(Surv(time, status) ~ A2BP1 + ABI3 + ACTL6B + ACVRL1 + ADAM11 + ADAM28 + ADAP1 + ADAP2 + ADORA3 + ADSL + AGAP4 + AGAP6 + AIF1 + AK5 + ALOX5 + ANG + ANKRD20A3 + ANKRD20B + ANKRD56 + ANXA2P1 + ANXA2P2 + ANXA2 + AOAH + APBB1IP + APLP1 + ARHGAP11A + ARHGAP30 + ARHGAP4 + ARHGAP9 + ARHGDIB + ARHGEF15 + ARL11 + ARL5B + ARPC1B + ASF1B + ASPM + ATAD3B + ATF2 + ATP1A3 + ATP5J2 + ATP5J + ATP5O + ATP8B4 + AURKA + AURKB + B3GNT4 + BCAS1 + BCL6B + BIN2 + BIRC5 + BIRC6 + BLNK + BOK + BRIP1 + BSN + BTF3 + BTK + BUB1B + BUB1 + C11orf87 + C11orf9 + C12orf59 + C16orf54 + C17orf60 + C19orf56 + C19orf70 + C1QA + C1QB + C1QC + C1QL3 + C1R + C1S + C1orf151 + C1orf173 + C1orf201 + C1orf38 + C20orf177 + C2orf49 + C3AR1 + C3 + C7orf59 + C9orf102 + CA7 + CACNG3 + CALY + CAMK1G + CAMK2A + CAPN3 + CARD9 + CARNS1 + CASC5 + CCDC75 + CCK + CCL3L1 + CCL3 + CCL4L2 + CCL4 + CCL5 + CCNA2 + CCNB2 + CCNT1 + CCR1 + CCR5 + CD14 + CD22 + CD2 + CD300A + CD300C + CD33 + CD37 + CD3E + CD4 + CD53 + CD68 + CD74 + CD86 + CDC20 + CDC25C +
# #                                       CDC45 + CDC6 + CDCA2 + CDCA5 + CDCA8 + CDK18 + CDK1 + CDK5R2 + CELF3 + CELF4 + CELF5 + CENPA + CENPE + CENPF + CENPK + CEP55 + CHD5 + CHGA + CHRM1 + CIITA + CKAP2L + CKMT1B + CLCA4 + CLDN11 + CLDND1 + CLEC18A + CLEC18B + CLEC2L + CLEC7A + CLIC1 + CLOCK + CLSPN + CMPK2 + CMTM7 + CNDP1 + CNP + CNTN2 + CNTNAP4 + COL1A1 + COL1A2 + COL3A1 + COL4A1 + COL4A2 + COL4A3 + COL4A4 + COL6A3 + CORO1A + COX6B1 + CPLX2 + CPNE6 + CPNE9 + CREG2 + CRYM + CSDAP1 + CSDA + CSF1R + CSF2RA + CSF2RB + CSF3R + CSNK2A1P + CSNK2A1 + CSRNP3 + CTSS + CXCR6 + CYB5R2 + CYBA + CYBB + CYTH4 + DAPP1 + DBNDD2 + DDI2 + DDN + DDX12 + DEF6 + DENND1C + DEPDC1B + DEPDC1 + DHRS4L2 + DHRS4 + DLG4 + DLGAP3 + DLGAP5 + DLX1 + DLX2 + DNAH17 + DNM1 + DOC2A + DOCK2 + DOCK8 + DTL + DTX3L + E2F2 + ECT2 + EEF1A1P9 + EEF1A1 + EEF1D + EEF1G + EIF3E + EIF3G + EIF3H + EIF3L + EIF5AL1 + EIF5A + ELOVL1 + EME1 + EMX1 + EMX2OS + EMX2 + ENPP2 + ENPP6 + EPHB6 + EPN1 + EPR1 + EPSTI1 + ERCC6L + ERMN + ESAM + ESCO2 + ESPL1 + EVI2B + EXO1 + EXOC6B + FA2H +
# #                                       FAM111B + FAM115C + FAM123C + FAM199X + FAM19A4 + FAM50B + FAM64A + FAM72B + FAM72D + FANCA + FANCD2 + FANCI + FAU + FBL + FBXL16 + FBXO41 + FCER1G + FCGR1A + FCGR1B + FCGR1C + FCGR2A + FCGR3A + FERMT3 + FER + FLT4 + FMNL1 + FOLH1 + FOS + FOXM1 + FOXO3B + FOXO3 + FPR1 + FYB + GABRA1 + GABRA5 + GABRB2 + GABRG2 + GAD2 + GDA + GJB1 + GJC2 + GLT1D1 + GMFG + GNA15 + GNB2L1 + GNG3 + GNGT2 + GPR22 + GPR26 + GPR34 + GPR37 + GPR62 + GPR65 + GPSM3 + GRAP + GRIN1 + GSG2 + GTF2A1 + GTF2IRD2P1 + GTF2IRD2 + GTF3C4 + GTSE1 + HAPLN2 + HAVCR2 + HBA2 + HBB + HCK + HCLS1 + HCN1 + HEPACAM + HEPN1 + HIPK4 + HJURP + HLAaaaA + HLAaaaDMA + HLAaaaDMB + HLAaaaDOA + HLAaaaDPA1 + HLAaaaDPB1 + HLAaaaDQA1 + HLAaaaDQB1 + HLAaaaDRA + HLAaaaDRB1 + HLAaaaDRB5 + HLAaaaH + HMHA1 + HMMR + HMP19 + HNRNPA1L2 + HNRNPA1 + HOOK1 + HOXD1 + HPCA + HPN + HRNBP3 + HSPA1A + HSPA1B + HTR5A + IFI30 + IFI6 + IFIT1 + IFIT3 + IKZF1 + IL10RA + IL12RB1 + IL18 + INA + IPMK + IQGAP3 + IQSEC3 + ISG15 + ISLR2 + ITGAL + ITGAM + ITGB2 + JUNB + KANK4 +
# #                                       KCNAB2 + KCNC2 + KCNK9 + KCNQ1DN + KCNS1 + KCNS2 + KCNT1 + KCNV1 + KIAA0101 + KIAA0748 + KIAA1045 + KIF11 + KIF14 + KIF20A + KIF23 + KIF2C + KIF4A + KIF4B + KIFC1 + KLHL26 + KLK6 + KPNA5 + LAIR1 + LAPTM5 + LAT2 + LATS1 + LCK + LCP1 + LDB3 + LDLRAP1 + LGALS9 + LGI3 + LHX6 + LILRB1 + LILRB4 + LMBRD2 + LMTK2 + LNPEP + LOC100132287 + LOC100133331 + LOC100233209 + LOC151009 + LOC154761 + LOC283999 + LOC285780 + LOC286002 + LOC441089 + LOC541471 + LOC606724 + LOC642846 + LOC653566 + LOC654433 + LPAR1 + LPAR5 + LRP2 + LRRC25 + LST1 + LUC7L + LY86 + LYL1 + MAG + MAL2 + MAL + MAP7D2 + MAPK8IP2 + MBP + MCM10 + MED27 + MELK + MKI67 + MLF1IP + MOBP + MOG + MPPED1 + MRAP2 + MRPS15 + MS4A4A + MS4A6A + MS4A7 + MSR1 + MX1 + MYBL2 + MYH6 + MYH7 + MYL12A + MYO1D + MYO1F + MYO9A + MYST3 + MYT1L + NACA + NAPB + NCAPG + NCAPH + NCDN + NCF4 + NCKAP1L + NCOA2 + NCOR1 + NCRNA00152 + NCRNA00188 + NDC80 + NDUFA13 + NDUFA2 + NDUFB7 + NECAB1 + NEFL + NEFM + NEIL3 + NEK2 + NEURL + NEUROD2 + NEUROD6 + NFKBIL2 + NHLRC2 +
# #                                       NINJ2 + NIPAL4 + NKAIN2 + NKX6aaa2 + NNAT + NRGN + NRIP3 + NSUN5P1 + NSUN5P2 + NUDT16P1 + NUF2 + NUSAP1 + OAS2 + OLFM3 + OLIG1 + OLIG2 + OLR1 + ORC1L + ORC6L + OSCAR + PABPC1L2A + PABPC1L2B + PABPC1 + PABPC3 + PACSIN1 + PANK3 + PARP9 + PARVG + PAX8 + PBK + PCP4L1 + PDIA3P + PDIA3 + PHF5A + PHYHIP + PIK3AP1 + PIK3CG + PIK3R5 + PILRB + PKMYT1 + PKN1 + PLCB2 + PLCG2 + PLD4 + PLEKHM3 + PLEK + PLK1 + PLK4 + PLP1 + PNMAL2 + POC1A + POLQ + POLR2I + POM121C + POM121 + PPAP2C + PPP1R14A + PRC1 + PREPL + PRKCG + PSD + PSMB8 + PSMB9 + PTAFR + PTER + PTGS1 + PTPN5 + PTPN6 + PTPRC + PTPRR + PYCARD + PYGO1 + RAB3A + RAB40B + RAD51 + RAD54L + RAPGEF5 + RASAL1 + RASAL3 + RASGRP4 + RBM47 + RBP4 + RBX1 + RC3H2 + RECQL4 + RGS10 + RGS18 + RGS19 + RGS4 + RGS7BP + RHOT2 + RIF1 + RIPK3 + RNASE2 + RNASE3 + RNASE4 + RNASE6 + RNASET2 + RNF11 + RPL10 + RPL11 + RPL12 + RPL13A + RPL13 + RPL14 + RPL15 + RPL18A + RPL18 + RPL19 + RPL23A + RPL23 + RPL24 + RPL27A + RPL27 + RPL29 + RPL30 + RPL31 + RPL32 + RPL34 + RPL35A + RPL36 +
# #                                       RPL37A + RPL37 + RPL3 + RPL4 + RPL5 + RPL6 + RPL7A + RPL7 + RPL8 + RPLP0 + RPLP2 + RPS10 + RPS11 + RPS12 + RPS13 + RPS15A + RPS15 + RPS16 + RPS18 + RPS19 + RPS20 + RPS21 + RPS23 + RPS24 + RPS25 + RPS27A + RPS3 + RPS4X + RPS5 + RPS6KA1 + RPS6 + RPS7 + RPS8 + RPSAP58 + RPSA + RRM2 + RSAD2 + RSC1A1 + RUNDC3A + RYR2 + S100A11 + S100A8 + S100A9 + S1PR5 + SAMD8 + SAMSN1 + SASH3 + SBNO1 + SCAF1 + SCIN + SCN2A + SCN8A + SCRT1 + SDR16C5 + SEC14L5 + SELPLG + SERPINI1 + SFRS16 + SGOL1 + SGOL2 + SH2D5 + SH3TC2 + SHCBP1 + SIGLEC7 + SIGLEC9 + SIPA1 + SKA1 + SKA3 + SLA + SLC12A5 + SLC15A3 + SLC17A7 + SLC2A14 + SLC2A3 + SLC2A5 + SLC30A3 + SLC31A2 + SLC32A1 + SLC45A3 + SLC5A11 + SLC6A17 + SLC6A7 + SLC7A7 + SLC8A2 + SMG1 + SNAP25 + SNAP91 + SNCB + SNRNP70 + SOCS7 + SPAG5 + SPC24 + SPC25 + SPCS2 + SPI1 + SPOCK3 + SRCIN1 + SRGN + SRRM4 + ST18 + STX1A + STX1B + STXBP1 + STXBP2 + SULT4A1 + SUMO1P3 + SUMO1 + SUSD3 + SV2B + SVOP + SYK + SYN1 + SYN2 + SYNGR2 + SYNGR3 + SYNJ1 + SYNJ2 + SYNPR + SYP + SYT13 + SYT1 + SYT5 +
# #                                       SYT7 + TACC3 + TAOK1 + TBCB + TBR1 + TBXAS1 + TFEC + TGFBRAP1 + TK1 + TLR1 + TLR2 + TLR7 + TM6SF1 + TMC8 + TMEM106A + TMEM119 + TMEM125 + TMEM130 + TMEM144 + TMEM151A + TMEM155 + TMEM176A + TMEM176B + TNFAIP8L2 + TOMM22 + TOP2A + TPX2 + TREM2 + TRHDE + TRIM59 + TROAP + TRPM2 + TTBK2 + TTC30A + TTC30B + TTK + TUBA8 + TUBB2A + TUBB2B + TUBB4 + TYROBP + UBE2C + UBE2MP1 + UBE2M + UBL5 + UBXN7 + UHMK1 + UNC13A + UNC13C + UQCR11 + UQCRHL + UQCRH + UQCRQ + USP34 + VAMP8 + VAV1 + VIP + VSIG4 + VSNL1 + WAS + WDFY4 + WDR47 + WDR62 + WNT10B + YIPF6 + ZCCHC12 + ZFP36 + ZFP57 + ZIC2 + ZIC5 + ZKSCAN1 + ZNF192 + ZNF354C + ZNF417 + ZNF426 + ZNF585A + ZNF585B + ZNF587 + ZNF678 + ZNF692 + ZNF699 + ZNF767 + ZNF778 + ZNF791 + ZNHIT1 + ZWINT,
# #                                     iter.max = 750, data = Pan_glioma_f_2016_and_SURV)
# # 
# # save(glasso_cox_model_643v_2016_it750, file = "~/paraRoberta/glasso_cox_model_643v_2016_it750.RData")
# # load(file = "~/paraRoberta/glasso_cox_model_643v_2016.RData")
# # 
# # #glasso_cox_model <- coxph(Surv(time, status) ~ A2BP1 + ABI3 + ACTL6B + ACVRL1 + ADAM11 + ADAM28 + ADAP1 + ADAP2 + ADORA3 + ADSL + AGAP4 + AGAP6 + AIF1 + AK5 + ALOX5 + ANG + ANKRD20A3 + ANKRD20B + ANKRD56 + ANXA2P1 + ANXA2P2 + ANXA2 + AOAH + APBB1IP + APLP1 + ARHGAP11A + ARHGAP30 + ARHGAP4 + ARHGAP9 + ARHGDIB + ARHGEF15 + ARL11 + ARL5B + ARPC1B + ASF1B + ASPM + ATAD3B + ATF2 + ATP1A3 + ATP5J2 + ATP5J + ATP5O + ATP8B4 + AURKA + AURKB + B3GNT4 + BCAS1 + BCL6B + BIN2 + BIRC5 + BIRC6 + BLNK + BOK + BRIP1 + BSN + BTF3 + BTK + BUB1B + BUB1 + C11orf87 + C11orf9 + C12orf59 + C16orf54 + C17orf60 + C19orf56 + C19orf70 + C1QA + C1QB + C1QC + C1QL3 + C1R + C1S + C1orf151 + C1orf173 + C1orf201 + C1orf38 + C20orf177 + C2orf49 + C3AR1 + C3 + C7orf59 + C9orf102 + CA7 + CACNG3 + CALY + CAMK1G + CAMK2A + CAPN3 + CARD9 + CARNS1 + CASC5 + CCDC75 + CCK + CCL3L1 + CCL3 + CCL4L2 + CCL4 + CCL5 + CCNA2 + CCNB2 + CCNT1 + CCR1 + CCR5 + CD14 + CD22 + CD2 + CD300A + CD300C + CD33 + CD37 + CD3E + CD4 + CD53 + CD68 + CD74 + CD86 + CDC20 + CDC25C + CDC45 + CDC6 + CDCA2 + CDCA5 + CDCA8 + CDK18 + CDK1 + CDK5R2 + CELF3 + CELF4 + CELF5 + CENPA + CENPE + CENPF + CENPK + CEP55 + CHD5 + CHGA + CHRM1 + CIITA + CKAP2L + CKMT1B + CLCA4 + CLDN11 + CLDND1 + CLEC18A + CLEC18B + CLEC2L + CLEC7A + CLIC1 + CLOCK + CLSPN + CMPK2 + CMTM7 + CNDP1 + CNP + CNTN2 + CNTNAP4 + COL1A1 + COL1A2 + COL3A1 + COL4A1 + COL4A2 + COL4A3 + COL4A4 + COL6A3 + CORO1A + COX6B1 + CPLX2 + CPNE6 + CPNE9 + CREG2 + CRYM + CSDAP1 + CSDA + CSF1R + CSF2RA + CSF2RB + CSF3R + CSNK2A1P + CSNK2A1 + CSRNP3 + CTSS + CXCR6 + CYB5R2 + CYBA + CYBB + CYTH4 + DAPP1 + DBNDD2 + DDI2 + DDN + DDX12 + DEF6 + DENND1C + DEPDC1B + DEPDC1 + DHRS4L2 + DHRS4 + DLG4 + DLGAP3 + DLGAP5 + DLX1 + DLX2 + DNAH17 + DNM1 + DOC2A + DOCK2 + DOCK8 + DTL + DTX3L + E2F2 + ECT2 + EEF1A1P9 + EEF1A1 + EEF1D + EEF1G + EIF3E + EIF3G + EIF3H + EIF3L + EIF5AL1 + EIF5A + ELOVL1 + EME1 + EMX1 + EMX2OS + EMX2 + ENPP2 + ENPP6 + EPHB6 + EPN1 + EPR1 + EPSTI1 + ERCC6L + ERMN + ESAM + ESCO2 + ESPL1 + EVI2B + EXO1 + EXOC6B + FA2H + FAM111B + FAM115C + FAM123C + FAM199X + FAM19A4 + FAM50B + FAM64A + FAM72B + FAM72D + FANCA + FANCD2 + FANCI + FAU + FBL + FBXL16 + FBXO41 + FCER1G + FCGR1A + FCGR1B + FCGR1C + FCGR2A + FCGR3A + FERMT3 + FER + FLT4 + FMNL1 + FOLH1 + FOS + FOXM1 + FOXO3B + FOXO3 + FPR1 + FYB + GABRA1 + GABRA5 + GABRB2 + GABRG2 + GAD2 + GDA + GJB1 + GJC2 + GLT1D1 + GMFG + GNA15 + GNB2L1 + GNG3 + GNGT2 + GPR22 + GPR26 + GPR34 + GPR37 + GPR62 + GPR65 + GPSM3 + GRAP + GRIN1 + GSG2 + GTF2A1 + GTF2IRD2P1 + GTF2IRD2 + GTF3C4 + GTSE1 + HAPLN2 + HAVCR2 + HBA2 + HBB + HCK + HCLS1 + HCN1 + HEPACAM + HEPN1 + HIPK4 + HJURP + HLA-A + HLA-DMA + HLA-DMB + HLA-DOA + HLA-DPA1 + HLA-DPB1 + HLA-DQA1 + HLA-DQB1 + HLA-DRA + HLA-DRB1 + HLA-DRB5 + HLA-H + HMHA1 + HMMR + HMP19 + HNRNPA1L2 + HNRNPA1 + HOOK1 + HOXD1 + HPCA + HPN + HRNBP3 + HSPA1A + HSPA1B + HTR5A + IFI30 + IFI6 + IFIT1 + IFIT3 + IKZF1 + IL10RA + IL12RB1 + IL18 + INA + IPMK + IQGAP3 + IQSEC3 + ISG15 + ISLR2 + ITGAL + ITGAM + ITGB2 + JUNB + KANK4 + KCNAB2 + KCNC2 + KCNK9 + KCNQ1DN + KCNS1 + KCNS2 + KCNT1 + KCNV1 + KIAA0101 + KIAA0748 + KIAA1045 + KIF11 + KIF14 + KIF20A + KIF23 + KIF2C + KIF4A + KIF4B + KIFC1 + KLHL26 + KLK6 + KPNA5 + LAIR1 + LAPTM5 + LAT2 + LATS1 + LCK + LCP1 + LDB3 + LDLRAP1 + LGALS9 + LGI3 + LHX6 + LILRB1 + LILRB4 + LMBRD2 + LMTK2 + LNPEP + LOC100132287 + LOC100133331 + LOC100233209 + LOC151009 + LOC154761 + LOC283999 + LOC285780 + LOC286002 + LOC441089 + LOC541471 + LOC606724 + LOC642846 + LOC653566 + LOC654433 + LPAR1 + LPAR5 + LRP2 + LRRC25 + LST1 + LUC7L + LY86 + LYL1 + MAG + MAL2 + MAL + MAP7D2 + MAPK8IP2 + MBP + MCM10 + MED27 + MELK + MKI67 + MLF1IP + MOBP + MOG + MPPED1 + MRAP2 + MRPS15 + MS4A4A + MS4A6A + MS4A7 + MSR1 + MX1 + MYBL2 + MYH6 + MYH7 + MYL12A + MYO1D + MYO1F + MYO9A + MYST3 + MYT1L + NACA + NAPB + NCAPG + NCAPH + NCDN + NCF4 + NCKAP1L + NCOA2 + NCOR1 + NCRNA00152 + NCRNA00188 + NDC80 + NDUFA13 + NDUFA2 + NDUFB7 + NECAB1 + NEFL + NEFM + NEIL3 + NEK2 + NEURL + NEUROD2 + NEUROD6 + NFKBIL2 + NHLRC2 + NINJ2 + NIPAL4 + NKAIN2 + NKX6-2 + NNAT + NRGN + NRIP3 + NSUN5P1 + NSUN5P2 + NUDT16P1 + NUF2 + NUSAP1 + OAS2 + OLFM3 + OLIG1 + OLIG2 + OLR1 + ORC1L + ORC6L + OSCAR + PABPC1L2A + PABPC1L2B + PABPC1 + PABPC3 + PACSIN1 + PANK3 + PARP9 + PARVG + PAX8 + PBK + PCP4L1 + PDIA3P + PDIA3 + PHF5A + PHYHIP + PIK3AP1 + PIK3CG + PIK3R5 + PILRB + PKMYT1 + PKN1 + PLCB2 + PLCG2 + PLD4 + PLEKHM3 + PLEK + PLK1 + PLK4 + PLP1 + PNMAL2 + POC1A + POLQ + POLR2I + POM121C + POM121 + PPAP2C + PPP1R14A + PRC1 + PREPL + PRKCG + PSD + PSMB8 + PSMB9 + PTAFR + PTER + PTGS1 + PTPN5 + PTPN6 + PTPRC + PTPRR + PYCARD + PYGO1 + RAB3A + RAB40B + RAD51 + RAD54L + RAPGEF5 + RASAL1 + RASAL3 + RASGRP4 + RBM47 + RBP4 + RBX1 + RC3H2 + RECQL4 + RGS10 + RGS18 + RGS19 + RGS4 + RGS7BP + RHOT2 + RIF1 + RIPK3 + RNASE2 + RNASE3 + RNASE4 + RNASE6 + RNASET2 + RNF11 + RPL10 + RPL11 + RPL12 + RPL13A + RPL13 + RPL14 + RPL15 + RPL18A + RPL18 + RPL19 + RPL23A + RPL23 + RPL24 + RPL27A + RPL27 + RPL29 + RPL30 + RPL31 + RPL32 + RPL34 + RPL35A + RPL36 + RPL37A + RPL37 + RPL3 + RPL4 + RPL5 + RPL6 + RPL7A + RPL7 + RPL8 + RPLP0 + RPLP2 + RPS10 + RPS11 + RPS12 + RPS13 + RPS15A + RPS15 + RPS16 + RPS18 + RPS19 + RPS20 + RPS21 + RPS23 + RPS24 + RPS25 + RPS27A + RPS3 + RPS4X + RPS5 + RPS6KA1 + RPS6 + RPS7 + RPS8 + RPSAP58 + RPSA + RRM2 + RSAD2 + RSC1A1 + RUNDC3A + RYR2 + S100A11 + S100A8 + S100A9 + S1PR5 + SAMD8 + SAMSN1 + SASH3 + SBNO1 + SCAF1 + SCIN + SCN2A + SCN8A + SCRT1 + SDR16C5 + SEC14L5 + SELPLG + SERPINI1 + SFRS16 + SGOL1 + SGOL2 + SH2D5 + SH3TC2 + SHCBP1 + SIGLEC7 + SIGLEC9 + SIPA1 + SKA1 + SKA3 + SLA + SLC12A5 + SLC15A3 + SLC17A7 + SLC2A14 + SLC2A3 + SLC2A5 + SLC30A3 + SLC31A2 + SLC32A1 + SLC45A3 + SLC5A11 + SLC6A17 + SLC6A7 + SLC7A7 + SLC8A2 + SMG1 + SNAP25 + SNAP91 + SNCB + SNRNP70 + SOCS7 + SPAG5 + SPC24 + SPC25 + SPCS2 + SPI1 + SPOCK3 + SRCIN1 + SRGN + SRRM4 + ST18 + STX1A + STX1B + STXBP1 + STXBP2 + SULT4A1 + SUMO1P3 + SUMO1 + SUSD3 + SV2B + SVOP + SYK + SYN1 + SYN2 + SYNGR2 + SYNGR3 + SYNJ1 + SYNJ2 + SYNPR + SYP + SYT13 + SYT1 + SYT5 + SYT7 + TACC3 + TAOK1 + TBCB + TBR1 + TBXAS1 + TFEC + TGFBRAP1 + TK1 + TLR1 + TLR2 + TLR7 + TM6SF1 + TMC8 + TMEM106A + TMEM119 + TMEM125 + TMEM130 + TMEM144 + TMEM151A + TMEM155 + TMEM176A + TMEM176B + TNFAIP8L2 + TOMM22 + TOP2A + TPX2 + TREM2 + TRHDE + TRIM59 + TROAP + TRPM2 + TTBK2 + TTC30A + TTC30B + TTK + TUBA8 + TUBB2A + TUBB2B + TUBB4 + TYROBP + UBE2C + UBE2MP1 + UBE2M + UBL5 + UBXN7 + UHMK1 + UNC13A + UNC13C + UQCR11 + UQCRHL + UQCRH + UQCRQ + USP34 + VAMP8 + VAV1 + VIP + VSIG4 + VSNL1 + WAS + WDFY4 + WDR47 + WDR62 + WNT10B + YIPF6 + ZCCHC12 + ZFP36 + ZFP57 + ZIC2 + ZIC5 + ZKSCAN1 + ZNF192 + ZNF354C + ZNF417 + ZNF426 + ZNF585A + ZNF585B + ZNF587 + ZNF678 + ZNF692 + ZNF699 + ZNF767 + ZNF778 + ZNF791 + ZNHIT1 + ZWINT, iter.max = 500, data =  Pan_glioma_f_2016_and_SURV)
# # #glasso_cox_model <- coxph(Surv(time, status) ~ A2BP1 + ABI3 + ACTL6B + ACVRL1 + ADAM11 + ADAM28 + ADAP1 + ADAP2 + ADORA3 + ADSL + AGAP4 + AGAP6 + AIF1 + AK5 + ALOX5 + ANG + ANKRD20A3 + ANKRD20B + ANKRD56 + ANXA2P1 + ANXA2P2 + ANXA2 + AOAH + APBB1IP + APLP1 + ARHGAP11A + ARHGAP30 + ARHGAP4 + ARHGAP9 + ARHGDIB + ARHGEF15 + ARL11 + ARL5B + ARPC1B + ASF1B + ASPM + ATAD3B + ATF2 + ATP1A3 + ATP5J2 + ATP5J + ATP5O + ATP8B4 + AURKA + AURKB + B3GNT4 + BCAS1 + BCL6B + BIN2 + BIRC5 + BIRC6 + BLNK + BOK + BRIP1 + BSN + BTF3 + BTK + BUB1B + BUB1 + C11orf87 + C11orf9 + C12orf59 + C16orf54 + C17orf60 + C19orf56 + C19orf70 + C1QA + C1QB + C1QC + C1QL3 + C1R + C1S + C1orf151 + C1orf173 + C1orf201 + C1orf38 + C20orf177 + C2orf49 + C3AR1 + C3 + C7orf59 + C9orf102 + CA7 + CACNG3 + CALY + CAMK1G + CAMK2A + CAPN3 + CARD9 + CARNS1 + CASC5 + CCDC75 + CCK + CCL3L1 + CCL3 + CCL4L2 + CCL4 + CCL5 + CCNA2 + CCNB2 + CCNT1 + CCR1 + CCR5 + CD14 + CD22 + CD2 + CD300A + CD300C + CD33 + CD37 + CD3E + CD4 + CD53 + CD68 + CD74 + CD86 + CDC20 + CDC25C + CDC45 + CDC6 + CDCA2 + CDCA5 + CDCA8 + CDK18 + CDK1 + CDK5R2 + CELF3 + CELF4 + CELF5 + CENPA + CENPE + CENPF + CENPK + CEP55 + CHD5 + CHGA + CHRM1 + CIITA + CKAP2L + CKMT1B + CLCA4 + CLDN11 + CLDND1 + CLEC18A + CLEC18B + CLEC2L + CLEC7A + CLIC1 + CLOCK + CLSPN + CMPK2 + CMTM7 + CNDP1 + CNP + CNTN2 + CNTNAP4 + COL1A1 + COL1A2 + COL3A1 + COL4A1 + COL4A2 + COL4A3 + COL4A4 + COL6A3 + CORO1A + COX6B1 + CPLX2 + CPNE6 + CPNE9 + CREG2 + CRYM + CSDAP1 + CSDA + CSF1R + CSF2RA + CSF2RB + CSF3R + CSNK2A1P + CSNK2A1 + CSRNP3 + CTSS + CXCR6 + CYB5R2 + CYBA + CYBB + CYTH4 + DAPP1 + DBNDD2 + DDI2 + DDN + DDX12 + DEF6 + DENND1C + DEPDC1B + DEPDC1 + DHRS4L2 + DHRS4 + DLG4 + DLGAP3 + DLGAP5 + DLX1 + DLX2 + DNAH17 + DNM1 + DOC2A + DOCK2 + DOCK8 + DTL + DTX3L + E2F2 + ECT2 + EEF1A1P9 + EEF1A1 + EEF1D + EEF1G + EIF3E + EIF3G + EIF3H + EIF3L + EIF5AL1 + EIF5A + ELOVL1 + EME1 + EMX1 + EMX2OS + EMX2 + ENPP2 + ENPP6 + EPHB6 + EPN1 + EPR1 + EPSTI1 + ERCC6L + ERMN + ESAM + ESCO2 + ESPL1 + EVI2B + EXO1 + EXOC6B + FA2H + FAM111B + FAM115C + FAM123C + FAM199X + FAM19A4 + FAM50B + FAM64A + FAM72B + FAM72D + FANCA + FANCD2 + FANCI + FAU + FBL + FBXL16 + FBXO41 + FCER1G + FCGR1A + FCGR1B + FCGR1C + FCGR2A + FCGR3A + FERMT3 + FER + FLT4 + FMNL1 + FOLH1 + FOS + FOXM1 + FOXO3B + FOXO3 + FPR1 + FYB + GABRA1 + GABRA5 + GABRB2 + GABRG2 + GAD2 + GDA + GJB1 + GJC2 + GLT1D1 + GMFG + GNA15 + GNB2L1 + GNG3 + GNGT2 + GPR22 + GPR26 + GPR34 + GPR37 + GPR62 + GPR65 + GPSM3 + GRAP + GRIN1 + GSG2 + GTF2A1 + GTF2IRD2P1 + GTF2IRD2 + GTF3C4 + GTSE1 + HAPLN2 + HAVCR2 + HBA2 + HBB + HCK + HCLS1 + HCN1 + HEPACAM + HEPN1 + HIPK4 + HJURP + HLA-A + HLA-DMA + HLA-DMB + HLA-DOA + HLA-DPA1 + HLA-DPB1 + HLA-DQA1 + HLA-DQB1 + HLA-DRA + HLA-DRB1 + HLA-DRB5 + HLA-H + HMHA1 + HMMR + HMP19 + HNRNPA1L2 + HNRNPA1 + HOOK1 + HOXD1 + HPCA + HPN + HRNBP3 + HSPA1A + HSPA1B + HTR5A + IFI30 + IFI6 + IFIT1 + IFIT3 + IKZF1 + IL10RA + IL12RB1 + IL18 + INA + IPMK + IQGAP3 + IQSEC3 + ISG15 + ISLR2 + ITGAL + ITGAM + ITGB2 + JUNB + KANK4 + KCNAB2 + KCNC2 + KCNK9 + KCNQ1DN + KCNS1 + KCNS2 + KCNT1 + KCNV1 + KIAA0101 + KIAA0748 + KIAA1045 + KIF11 + KIF14 + KIF20A + KIF23 + KIF2C + KIF4A + KIF4B + KIFC1 + KLHL26 + KLK6 + KPNA5 + LAIR1 + LAPTM5 + LAT2 + LATS1 + LCK + LCP1 + LDB3 + LDLRAP1 + LGALS9 + LGI3 + LHX6 + LILRB1 + LILRB4 + LMBRD2 + LMTK2 + LNPEP + LOC100132287 + LOC100133331 + LOC100233209 + LOC151009 + LOC154761 + LOC283999 + LOC285780 + LOC286002 + LOC441089 + LOC541471 + LOC606724 + LOC642846 + LOC653566 + LOC654433 + LPAR1 + LPAR5 + LRP2 + LRRC25 + LST1 + LUC7L + LY86 + LYL1 + MAG + MAL2 + MAL + MAP7D2 + MAPK8IP2 + MBP + MCM10 + MED27 + MELK + MKI67 + MLF1IP + MOBP + MOG + MPPED1 + MRAP2 + MRPS15 + MS4A4A + MS4A6A + MS4A7 + MSR1 + MX1 + MYBL2 + MYH6 + MYH7 + MYL12A + MYO1D + MYO1F + MYO9A + MYST3 + MYT1L + NACA + NAPB + NCAPG + NCAPH + NCDN + NCF4 + NCKAP1L + NCOA2 + NCOR1 + NCRNA00152 + NCRNA00188 + NDC80 + NDUFA13 + NDUFA2 + NDUFB7 + NECAB1 + NEFL + NEFM + NEIL3 + NEK2 + NEURL + NEUROD2 + NEUROD6 + NFKBIL2 + NHLRC2 + NINJ2 + NIPAL4 + NKAIN2 + NKX6-2 + NNAT + NRGN + NRIP3 + NSUN5P1 + NSUN5P2 + NUDT16P1 + NUF2 + NUSAP1 + OAS2 + OLFM3 + OLIG1 + OLIG2 + OLR1 + ORC1L + ORC6L + OSCAR + PABPC1L2A + PABPC1L2B + PABPC1 + PABPC3 + PACSIN1 + PANK3 + PARP9 + PARVG + PAX8 + PBK + PCP4L1 + PDIA3P + PDIA3 + PHF5A + PHYHIP + PIK3AP1 + PIK3CG + PIK3R5 + PILRB + PKMYT1 + PKN1 + PLCB2 + PLCG2 + PLD4 + PLEKHM3 + PLEK + PLK1 + PLK4 + PLP1 + PNMAL2 + POC1A + POLQ + POLR2I + POM121C + POM121 + PPAP2C + PPP1R14A + PRC1 + PREPL + PRKCG + PSD + PSMB8 + PSMB9 + PTAFR + PTER + PTGS1 + PTPN5 + PTPN6 + PTPRC + PTPRR + PYCARD + PYGO1 + RAB3A + RAB40B + RAD51 + RAD54L + RAPGEF5 + RASAL1 + RASAL3 + RASGRP4 + RBM47 + RBP4 + RBX1 + RC3H2 + RECQL4 + RGS10 + RGS18 + RGS19 + RGS4 + RGS7BP + RHOT2 + RIF1 + RIPK3 + RNASE2 + RNASE3 + RNASE4 + RNASE6 + RNASET2 + RNF11 + RPL10 + RPL11 + RPL12 + RPL13A + RPL13 + RPL14 + RPL15 + RPL18A + RPL18 + RPL19 + RPL23A + RPL23 + RPL24 + RPL27A + RPL27 + RPL29 + RPL30 + RPL31 + RPL32 + RPL34 + RPL35A + RPL36 + RPL37A + RPL37 + RPL3 + RPL4 + RPL5 + RPL6 + RPL7A + RPL7 + RPL8 + RPLP0 + RPLP2 + RPS10 + RPS11 + RPS12 + RPS13 + RPS15A + RPS15 + RPS16 + RPS18 + RPS19 + RPS20 + RPS21 + RPS23 + RPS24 + RPS25 + RPS27A + RPS3 + RPS4X + RPS5 + RPS6KA1 + RPS6 + RPS7 + RPS8 + RPSAP58 + RPSA + RRM2 + RSAD2 + RSC1A1 + RUNDC3A + RYR2 + S100A11 + S100A8 + S100A9 + S1PR5 + SAMD8 + SAMSN1 + SASH3 + SBNO1 + SCAF1 + SCIN + SCN2A + SCN8A + SCRT1 + SDR16C5 + SEC14L5 + SELPLG + SERPINI1 + SFRS16 + SGOL1 + SGOL2 + SH2D5 + SH3TC2 + SHCBP1 + SIGLEC7 + SIGLEC9 + SIPA1 + SKA1 + SKA3 + SLA + SLC12A5 + SLC15A3 + SLC17A7 + SLC2A14 + SLC2A3 + SLC2A5 + SLC30A3 + SLC31A2 + SLC32A1 + SLC45A3 + SLC5A11 + SLC6A17 + SLC6A7 + SLC7A7 + SLC8A2 + SMG1 + SNAP25 + SNAP91 + SNCB + SNRNP70 + SOCS7 + SPAG5 + SPC24 + SPC25 + SPCS2 + SPI1 + SPOCK3 + SRCIN1 + SRGN + SRRM4 + ST18 + STX1A + STX1B + STXBP1 + STXBP2 + SULT4A1 + SUMO1P3 + SUMO1 + SUSD3 + SV2B + SVOP + SYK + SYN1 + SYN2 + SYNGR2 + SYNGR3 + SYNJ1 + SYNJ2 + SYNPR + SYP + SYT13 + SYT1 + SYT5 + SYT7 + TACC3 + TAOK1 + TBCB + TBR1 + TBXAS1 + TFEC + TGFBRAP1 + TK1 + TLR1 + TLR2 + TLR7 + TM6SF1 + TMC8 + TMEM106A + TMEM119 + TMEM125 + TMEM130 + TMEM144 + TMEM151A + TMEM155 + TMEM176A + TMEM176B + TNFAIP8L2 + TOMM22 + TOP2A + TPX2 + TREM2 + TRHDE + TRIM59 + TROAP + TRPM2 + TTBK2 + TTC30A + TTC30B + TTK + TUBA8 + TUBB2A + TUBB2B + TUBB4 + TYROBP + UBE2C + UBE2MP1 + UBE2M + UBL5 + UBXN7 + UHMK1 + UNC13A + UNC13C + UQCR11 + UQCRHL + UQCRH + UQCRQ + USP34 + VAMP8 + VAV1 + VIP + VSIG4 + VSNL1 + WAS + WDFY4 + WDR47 + WDR62 + WNT10B + YIPF6 + ZCCHC12 + ZFP36 + ZFP57 + ZIC2 + ZIC5 + ZKSCAN1 + ZNF192 + ZNF354C + ZNF417 + ZNF426 + ZNF585A + ZNF585B + ZNF587 + ZNF678 + ZNF692 + ZNF699 + ZNF767 + ZNF778 + ZNF791 + ZNHIT1 + ZWINT, iter.max = 500, data =  Pan_glioma_f_2016_and_SURV)
# # #glasso_cox_model <- coxph(Surv(time, status) ~ A2BP1 + PTPRR + PYCARD + PYGO1 + RAB3A + RAB40B + RAD51 + RAD54L + RAPGEF5 + RASAL1 + RASAL3 + RASGRP4 + RBM47 + RBP4 + RBX1 + RC3H2 + RECQL4 + RGS10 + RGS18 + RGS19 + RGS4 + RGS7BP + RHOT2 + RIF1 + RIPK3 + RNASE2 + RNASE3 + RNASE4 + RNASE6 + RNASET2 + RNF11 + RPL10 + RPL11 + RPL12 + RPL13A + RPL13 + RPL14 + RPL15 + RPL18A + RPL18 + RPL19 + RPL23A + RPL23 + RPL24 + RPL27A + RPL27 + RPL29 + RPL30 + RPL31 + RPL32 + RPL34 + RPL35A + RPL36 + RPL37A + RPL37 + RPL3 + RPL4 + RPL5 + RPL6 + RPL7A + RPL7 + RPL8 + RPLP0 + RPLP2 + RPS10 + RPS11 + RPS12 + RPS13 + RPS15A + RPS15 + RPS16 + RPS18 + RPS19 + RPS20 + RPS21 + RPS23 + RPS24 + RPS25 + RPS27A + RPS3 + RPS4X + RPS5 + RPS6KA1 + RPS6 + RPS7 + RPS8 + RPSAP58 + RPSA + RRM2 + RSAD2 + RSC1A1 + RUNDC3A + RYR2 + S100A11 + S100A8 + S100A9 + S1PR5 + SAMD8 + SAMSN1 + SASH3 + SBNO1 + SCAF1 + SCIN + SCN2A + SCN8A + SCRT1 + SDR16C5 + SEC14L5 + SELPLG + SERPINI1 + SFRS16 + SGOL1 + SGOL2 + SH2D5 + SH3TC2 + SHCBP1 + SIGLEC7 + SIGLEC9 + SIPA1 + SKA1 + SKA3 + SLA + SLC12A5 + SLC15A3 + SLC17A7 + SLC2A14 + SLC2A3 + SLC2A5 + SLC30A3 + SLC31A2 + SLC32A1 + SLC45A3 + SLC5A11 + SLC6A17 + SLC6A7 + SLC7A7 + SLC8A2 + SMG1 + SNAP25 + SNAP91 + SNCB + SNRNP70 + SOCS7 + SPAG5 + SPC24 + SPC25 + SPCS2 + SPI1 + SPOCK3 + SRCIN1 + SRGN + SRRM4 + ST18 + STX1A + STX1B + STXBP1 + STXBP2 + SULT4A1 + SUMO1P3 + SUMO1 + SUSD3 + SV2B + SVOP + SYK + SYN1 + SYN2 + SYNGR2 + SYNGR3 + SYNJ1 + SYNJ2 + SYNPR + SYP + SYT13 + SYT1 + SYT5 + SYT7 + TACC3 + TAOK1 + TBCB + TBR1 + TBXAS1 + TFEC + TGFBRAP1 + TK1 + TLR1 + TLR2 + TLR7 + TM6SF1 + TMC8 + TMEM106A + TMEM119 + TMEM125 + TMEM130 + TMEM144 + TMEM151A + TMEM155 + TMEM176A + TMEM176B + TNFAIP8L2 + TOMM22 + TOP2A + TPX2 + TREM2 + TRHDE + TRIM59 + TROAP + TRPM2 + TTBK2 + TTC30A + TTC30B + TTK + TUBA8 + TUBB2A + TUBB2B + TUBB4 + TYROBP + UBE2C + UBE2MP1 + UBE2M + UBL5 + UBXN7 + UHMK1 + UNC13A + UNC13C + UQCR11 + UQCRHL + UQCRH + UQCRQ + USP34 + VAMP8 + VAV1 + VIP + VSIG4 + VSNL1 + WAS + WDFY4 + WDR47 + WDR62 + WNT10B + YIPF6 + ZCCHC12 + ZFP36 + ZFP57 + ZIC2 + ZIC5 + ZKSCAN1 + ZNF192 + ZNF354C + ZNF417 + ZNF426 + ZNF585A + ZNF585B + ZNF587 + ZNF678 + ZNF692 + ZNF699 + ZNF767 + ZNF778 + ZNF791 + ZNHIT1 + ZWINT, iter.max = 500, data =  Pan_glioma_f_2016_and_SURV)
# # 
# # #library(survivalAnalysis)
# # summary_cox_it500     <- summary(glasso_cox_model_643v_2016_it500)
# # df_summary_cox_it500  <- as.data.frame(summary_cox_it500$coefficients)
# # #df_summary <- cox_as_data_frame(summary_cox_it500)
# # #summary(res.cox0000000001)
# 
# 
# 
# 
# 
# # A.2 - Cox model with features selected from GLASSO (2016: 643 variables) + LASSO (alpha=1)
# glasso_lasso_cox_model_643v_2016 <- glmnet(Panglioma_2016_ALL, Surv(SURV$time, SURV$status), alpha = 1, family = "cox")
# 
# #Surv(SURV$time, SURV$status)
# #sum(SURV$status==1) #EVENTS = 233
# #RULE = 15 / 10 EPV
# #pmax: Limit the maximum number of variables ever to be nonzero
# # events = 233 (deaths)
# #pmax = 233 / 15 = 15
# #pmax = 233 / 10 = 23
# 
# rm(Pan_glioma_i_2016_2, Panglioma_2016_ALL, Pan_glioma_i_2016_and_SURV,
#    Pan_glioma_f_2016_and_SURV, glasso_lasso_cox_model_643v_2016)



# 10 EPV ------------------------------------------------------------------
# = 233 / 10 = 23,3
Panglioma_2016_ALL <-Panglioma_2016_ALL[complete.cases(Panglioma_2016_ALL),]

glasso_lasso_cox_model_643v_2016_10EPV <- glmnet(Panglioma_2016_ALL, 
                                           Surv(SURV$time, SURV$status),
                                           alpha = 1, 
                                           family = "cox",
                                           pmax = 30) #pmax = 24 XXX

print(glasso_lasso_cox_model_643v_2016_10EPV)


ix_2016_10EPV           <- length(glasso_lasso_cox_model_643v_2016_10EPV$lambda)
lambda_value_2016_10EPV <- glasso_lasso_cox_model_643v_2016_10EPV$lambda[ix_2016_10EPV]
model_2016_10EPV        <- glmnet(Panglioma_2016_ALL, Surv(SURV$time, SURV$status),
                                  alpha = 1,
                                  lambda = lambda_value_2016_10EPV,
                                  family = "cox")

# extract coefficients at a single value of lambda
CF_2016_10EPV <- as.matrix(coef(glasso_lasso_cox_model_643v_2016_10EPV, s = lambda_value_2016_10EPV))  
CF_2016_10EPV <- as.data.frame(CF_2016_10EPV[CF_2016_10EPV!=0,])
CF_2016_10EPV <- CF_2016_10EPV[order(CF_2016_10EPV[1,])]
names(CF_2016_10EPV) <- "coef"
print(CF_2016_10EPV)

#string / formula for Kaplan Meier formula (afterwards)
#STRING_CF_2016_10EPV <- paste0(rownames(CF_2016_10EPV), collapse = " + ")
#STRING_CF_2016_10EPV

### MODELLING - GLMSPARSENET - SEPARATE2GROUPS
#USING COEFS FROM ELASTIC NET PERFORMED ON ALL 2016 INITIAL DATASET
CF_2016_10EPV <- as.matrix(coef(glasso_lasso_cox_model_643v_2016_10EPV, 
                     s = lambda_value_2016_10EPV))

sep2groups2016_10EPV <- glmSparseNet::separate2GroupsCox(as.vector(CF_2016_10EPV), 
                                                   Panglioma_2016_ALL[, rownames(CF_2016_10EPV)], 
                                                   SURV,
                                                   probs = c(0.5, 0.5),
                                                   plot.title = 'GLASSO+LASSO (2016_10EPV = 24 features)', legend.outside = FALSE)

#PLOTS, UNDERSTAND PROGNOSTIC INDEX
table(sep2groups2016_10EPV[["km"]][["custom.data"]][["group"]])
length(sep2groups2016_10EPV[["km"]][["custom.data"]][["group"]])


#############################################################
############### PROGNOSTIC INDEX: DISTRIBUTION ##############
#############################################################

# NEW - after meeting at NOVA FCT

#GRAPHICAL APPROACH
#BASED ON DENSITY FUNCTION ESTIMATION AND FINDING LOCAL MINIMUM (after meeting at NOVA)
#sep2groups2016_10EPV[["km"]][["custom.data"]][["group"]]
dat <- data.frame(sep2groups2016_10EPV$km$custom.data$group,sep2groups2016_10EPV$km$custom.data$pi)

# AFTER MCLUST, LOOK AGAIN AT P.I. distribution and
hist(sep2groups2016_10EPV$km$custom.data$pi)
GEOM_DENSITY_PLOT <- ggplot(dat, aes(x=sep2groups2016_10EPV$km$custom.data$pi)) + 
  geom_density() + 
  labs(title="2016, single omics, first look at PI to infer the no. of mixture components (clusters)")
#+ geom_vline(aes(xintercept=median(sep2groups2016_10EPV$km$custom.data$pi, na.rm=T)),color="red", linetype="dashed", size=1
GEOM_DENSITY_PLOT



library(mclust)
#Model-based clustering based on parameterized finite Gaussian mixture models.
#Models are estimated by EM algorithm initialized by hierarchical model-based agglomerative clustering.
#The optimal model is then selected according to BIC.

mclust_2 <- Mclust(sep2groups2016_10EPV$km$custom.data$pi, G = 2, modelNames = "V", 
                   prior = NULL, control = emControl(), initialization = NULL, 
                   warn = mclust.options("warn"), x =  NULL, verbose = interactive())
#output
mclust_2[["parameters"]][["mean"]]
mclust_2[["parameters"]][["mean"]][["1"]]
mclust_2[["parameters"]][["mean"]][["2"]]

mclust_3 <- Mclust(sep2groups2016_10EPV$km$custom.data$pi, G = 3, modelNames = "V", 
                   prior = NULL, control = emControl(), initialization = NULL, 
                   warn = mclust.options("warn"), x =  NULL, verbose = interactive())
#output
mclust_3[["parameters"]][["mean"]]
mclust_3[["parameters"]][["mean"]][["1"]]
mclust_3[["parameters"]][["mean"]][["2"]]
mclust_3[["parameters"]][["mean"]][["3"]]


mclust_4 <- Mclust(sep2groups2016_10EPV$km$custom.data$pi, G = 4, modelNames = "V", 
                   prior = NULL, control = emControl(), initialization = NULL, 
                   warn = mclust.options("warn"), x =  NULL, verbose = interactive())

#output
mclust_4[["parameters"]][["mean"]]
mclust_4[["parameters"]][["mean"]][["1"]]
mclust_4[["parameters"]][["mean"]][["2"]]
mclust_4[["parameters"]][["mean"]][["3"]]
mclust_4[["parameters"]][["mean"]][["4"]]

vertical.lines_mclust2 <- c(mclust_2[["parameters"]][["mean"]][["1"]],
                            mclust_2[["parameters"]][["mean"]][["2"]])

# AFTER MCLUST, LOOK AGAIN AT P.I. distribution and
GEOM_DENSITY_PLOT <- ggplot(dat, aes(x=sep2groups2016_10EPV$km$custom.data$pi)) + geom_density() 
#+ geom_vline(aes(xintercept=median(sep2groups2016_10EPV$km$custom.data$pi, na.rm=T)),color="red", linetype="dashed", size=1

# VERTICAL LINES FROM MCLUST_G=2 RESULTS
ggplot(dat, aes(x=sep2groups2016_10EPV$km$custom.data$pi)) + 
  labs(title="2016, single omics \noutput when no. of clusters = 2 \n(mean values of estimated distributions)") +
  geom_density() + geom_vline(xintercept = vertical.lines_mclust2, linetype="dotted", size = 0.3) +
  annotate("text", x = vertical.lines[1], y = 0.3, label = round(vertical.lines_mclust2[1],4), angle=90) +
  annotate("text", x = vertical.lines[2], y = 0.3, label = round(vertical.lines_mclust2[2],4), color="red", angle=90)
  



#MANUALLY CALCULATE PERCENTAGE FOR GROUP STRATIFICATION
# with the most suitable value from mclust output
pi_vector <- data.frame(sep2groups2016_10EPV$km$custom.data$pi)

#length(pi_vector[pi_vector<=0.397])
length(pi_vector[pi_vector<= mclust_2[["parameters"]][["mean"]][["2"]]]) #VALUE BASED ON MCLUST OUTPUT?
dim(pi_vector)[1]
percentil_2016_mclust2 <- length(pi_vector[pi_vector<= mclust_2[["parameters"]][["mean"]][["2"]]]) / dim(pi_vector)[1]


#percentil to stratify into low-/high-risk group
cat("INTERMEDIATE RESULT \n2016 single omics, Panglioma=LGG vs. GBM\nWe will divide high-/low-risk groups with percentil = ", percentil_2016_mclust2)

sep2groups2016_10EPV <- glmSparseNet::separate2GroupsCox(as.vector(CF_2016_10EPV), 
                                                         Panglioma_2016_ALL[, rownames(CF_2016_10EPV)], 
                                                         SURV,
                                                         probs = c(percentil_2016_mclust2, percentil_2016_mclust2),
                                                         plot.title = 'GLASSO+LASSO (2016_10EPV = 24 features) AFTER mclust', legend.outside = FALSE)
sep2groups2016_10EPV

density_df <- density(sep2groups2016_10EPV$km$custom.data$pi, n=2048)#n=1024)
density_df <- data.frame(density_df$x, density_df$y)


#library(spatialEco)
#local minimum
#abline(v=dx$x[localMinima(dx$y)],col=2,lty=2)
#lmm <- local.min.max(density_df$x, dev=mean, add.points=FALSE, 
#                     main="Local Minima and Maxima")
#plot(density_df)
#abline(v=density_df$density_df.x[localMinima(density_df$density_df.y)],col=2,lty=2)




#PLOT DECONVOLUTED GAUSSIANS (2)
# sigmasq_1 <- mclust[["parameters"]][["variance"]][["sigmasq"]][1]
# mean_1      <- mclust[["parameters"]][["mean"]][["1"]]
# sigmasq_2 <- mclust[["parameters"]][["variance"]][["sigmasq"]][2]
# mean_2      <- mclust[["parameters"]][["mean"]][["2"]]
# 
# #Create a sequence of 1000 x values based on population mean and standard deviation
# x_aux_1 <- seq(-3.1647, 3.1647, length = 1000) * sqrt(sigmasq_1) + mean_1
# y_aux_1 <- dnorm(x_aux_1, mean_1, sqrt(sigmasq_1))
# 
# x_aux_1 <- 
#   
# x_aux_1 <- rnorm(n= 1000, mean_1,sqrt(sigmasq_1))
# y_aux_1 <- dnorm(x_aux_1, mean_1, sqrt(sigmasq_1))
# 
# 
# #plot normal distribution with customized x-axis labels
# plot(x_aux_1,y_aux_1, type = "l", lwd = 2, axes = FALSE, xlab = "", ylab = "")
# sd_axis_bounds = 5
# axis_bounds <- seq(-sd_axis_bounds * sigmasq_1 + mean_1,
#                    sd_axis_bounds * sigmasq_1 + mean_1,
#                    by = sigmasq_1)
# axis(side = 1, at = axis_bounds, pos = 0)
# 
# plot(x_aux_1,y_aux_1, type = "l", lwd = 2, axes = TRUE, xlab = "", ylab = "")
# #axis(1, at = -3:3, labels = c("-3s", "-2s", "-1s", "mean", "1s", "2s", "3s"))
# 
# x_aux_1 <- seq(-4, 4, length = 1000) * sqrt(sigmasq_1) + mean_1
# y_aux_1 <- dnorm(x_aux_1, mean_1, sqrt(sigmasq_1))
# 
# 
# x <- sep2groups2016_10EPV$km$custom.data$pi
# hist(x)


#AUXILIARY
#STRING_CF_2016_10EPV
# fit <- survfit(Surv(time, status) ~ AGAP4 + AURKA + CDC6 + CLEC18A + CLEC18B + DHRS4 + EEF1A1P9 + EIF3L +
#                  FAM115C + GTF2IRD2 + HMP19 + KIF20A + MYBL2 + OSCAR + PABPC3 + POM121C + RPL7A + RRM2 +
#                  S100A11 + SAMD8 + SGOL1 + SNAP91 + TMEM176A + USP34, data = Pan_glioma_f_2016_and_SURV[, rownames(CF_2016_10EPV)])
# 
# ggsurvplot(fit,
#            pval = TRUE, conf.int = TRUE,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            linetype = "strata", # Change line type by groups
#            #surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            palette = c("#E7B800", "#2E9FDF"))
# 
# surv_diff <- survdiff(Surv(time, status) ~ radiotherapy, data = df_merged_analysis3)
# surv_diff


###################

#DIVIDE INTO LOW AND HIGH GROUPS
d_10EPV_2016             <- data.frame(sep2groups2016_10EPV[["km"]][["custom.data"]])
d_10EPV_2016$Patient_ID <- rownames(d_10EPV_2016)

d_10EPV_2016 <- merge(d_10EPV_2016,panglioma_clinic_single_omics_2016[,c(1,11)],by="Patient_ID")

#OUTPUT LIST
df_low_10EPV_2016  <- d_10EPV_2016[d_10EPV_2016$group == "Low risk", ]
df_high_10EPV_2016 <- d_10EPV_2016[d_10EPV_2016$group == "High risk", ]

write.csv(df_low_10EPV_2016,"~/paraRoberta/single-omics-ds/list_LOWgroup_10EPV_2016_mclust2.csv", row.names = TRUE)
write.csv(df_high_10EPV_2016,"~/paraRoberta/single-omics-ds/list_HIGHgroup_10EPV_2016_mclust2.csv", row.names = TRUE)

# STATISTICS / PERCENTAGES
#AUX
low_group_patients_10EPV_2016  <- df_low_10EPV_2016$Patient_ID
high_group_patients_10EPV_2016 <- df_high_10EPV_2016$Patient_ID

#table(panglioma_clinic_single_omics_2016$histological_type)

#####old
lgg_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type %in% c('astrocytoma','oligoastrocytoma','oligodendroglioma'),]
lgg_patients <- lgg_patients$Patient_ID

####new, after meeting
astro_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type %in% c('astrocytoma'),]
astro_patients <- astro_patients$Patient_ID

oligo_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type %in% c('oligodendroglioma'),]
oligo_patients <- oligo_patients$Patient_ID

mix_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type %in% c('oligoastrocytoma'),]
mix_patients <- mix_patients$Patient_ID


gbm_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type == "untreated primary (de novo) gbm",]
gbm_patients <- gbm_patients$Patient_ID

length(lgg_patients) #494 LGG
length(astro_patients) #193 LGG
length(oligo_patients) #188 LGG
length(mix_patients) #113 LGG
length(gbm_patients) #149 GBM


# LOW RISK GROUP
#GBM
percent_GBM_LOWrisk_group <- (sum(low_group_patients_10EPV_2016 %in% gbm_patients) * 100 ) / length(low_group_patients_10EPV_2016)
round(percent_GBM_LOWrisk_group,1)

#LGG
percent_LGG_LOWrisk_group <- (sum(low_group_patients_10EPV_2016 %in% lgg_patients) * 100 ) / length(low_group_patients_10EPV_2016)
round(percent_LGG_LOWrisk_group,1)

percent_ASTRO_LOWrisk_group <- (sum(low_group_patients_10EPV_2016 %in% astro_patients) * 100 ) / length(low_group_patients_10EPV_2016)
round(percent_ASTRO_LOWrisk_group,1)

percent_OLIGO_LOWrisk_group <- (sum(low_group_patients_10EPV_2016 %in% oligo_patients) * 100 ) / length(low_group_patients_10EPV_2016)
round(percent_OLIGO_LOWrisk_group,1)

percent_MIX_LOWrisk_group <- (sum(low_group_patients_10EPV_2016 %in% mix_patients) * 100 ) / length(low_group_patients_10EPV_2016)
round(percent_MIX_LOWrisk_group,1)

#check
percent_ASTRO_LOWrisk_group + percent_OLIGO_LOWrisk_group + percent_MIX_LOWrisk_group + percent_GBM_LOWrisk_group

##
# HIGH RISK GROUP
#GBM
percent_GBM_HIGHrisk_group <- (sum(high_group_patients_10EPV_2016 %in% gbm_patients) * 100 ) / length(high_group_patients_10EPV_2016)
round(percent_GBM_HIGHrisk_group,1)

#LGG
percent_LGG_HIGHrisk_group <- (sum(high_group_patients_10EPV_2016 %in% lgg_patients) * 100 ) / length(high_group_patients_10EPV_2016)
round(percent_LGG_HIGHrisk_group,1)

percent_ASTRO_HIGHrisk_group <- (sum(high_group_patients_10EPV_2016 %in% astro_patients) * 100 ) / length(high_group_patients_10EPV_2016)
round(percent_ASTRO_HIGHrisk_group,1)

percent_OLIGO_HIGHrisk_group <- (sum(high_group_patients_10EPV_2016 %in% oligo_patients) * 100 ) / length(high_group_patients_10EPV_2016)
round(percent_OLIGO_HIGHrisk_group,1)

percent_MIX_HIGHrisk_group <- (sum(high_group_patients_10EPV_2016 %in% mix_patients) * 100 ) / length(high_group_patients_10EPV_2016)
round(percent_MIX_HIGHrisk_group,1)

#check
percent_ASTRO_HIGHrisk_group + percent_OLIGO_HIGHrisk_group + percent_MIX_HIGHrisk_group + percent_GBM_HIGHrisk_group
percent_ASTRO_HIGHrisk_group + percent_OLIGO_HIGHrisk_group + percent_MIX_HIGHrisk_group

cat("% astrocytomas in high-risk group:", percent_ASTRO_HIGHrisk_group)
cat("% oligodendrogliomas in high-risk group:", percent_OLIGO_HIGHrisk_group)
cat("% mix (oligoastrocytomas) in high-risk group:", percent_MIX_HIGHrisk_group)
cat("% gbm in high-risk group:", percent_GBM_HIGHrisk_group)



#AFTER NOVA FCT MEETING
#low_group_patients_10EPV_2016 %in% gbm_patients
sum(low_group_patients_10EPV_2016 %in% gbm_patients)
intersect(low_group_patients_10EPV_2016, gbm_patients)
#OUTPUT LIST
lowriskgroup_patients_that_are_gbm_10EPV_2016 <- data.frame(d_10EPV_2016[d_10EPV_2016$Patient_ID %in% intersect(low_group_patients_10EPV_2016, gbm_patients), ])

sum(high_group_patients_10EPV_2016 %in% lgg_patients)
intersect(high_group_patients_10EPV_2016, lgg_patients)
#OUTPUT LIST
highriskgroup_patients_that_are_lgg_10EPV_2016 <- data.frame(d_10EPV_2016[d_10EPV_2016$Patient_ID %in% intersect(high_group_patients_10EPV_2016, lgg_patients), ])



write.csv(lowriskgroup_patients_that_are_gbm_10EPV_2016,"~/paraRoberta/single-omics-ds/list_LOWgroup_that_are_GBM_10EPV_2016_mclust2.csv", row.names = TRUE)
write.csv(highriskgroup_patients_that_are_lgg_10EPV_2016,"~/paraRoberta/single-omics-ds/list_HIGHgroup_that_are_LGG_10EPV_2016_mclust2.csv", row.names = TRUE)






























# 12 EPV ------------------------------------------------------------------
# = 233 / 12 = 19.41

glasso_lasso_cox_model_643v_2016_12EPV <- glmnet(Panglioma_2016_ALL, 
                                                 Surv(SURV$time, SURV$status),
                                                 alpha = 1, 
                                                 family = "cox",
                                                 pmax = 25) #pmax = 24 XXX

print(glasso_lasso_cox_model_643v_2016_12EPV)


ix_2016_12EPV           <- length(glasso_lasso_cox_model_643v_2016_12EPV$lambda)
lambda_value_2016_12EPV <- glasso_lasso_cox_model_643v_2016_12EPV$lambda[ix_2016_12EPV]
model_2016_12EPV        <- glmnet(Panglioma_2016_ALL, Surv(SURV$time, SURV$status),
                                  alpha = 1,
                                  lambda = lambda_value_2016_12EPV,
                                  family = "cox")

# extract coefficients at a single value of lambda
CF_2016_12EPV <- as.matrix(coef(glasso_lasso_cox_model_643v_2016_12EPV, s = lambda_value_2016_12EPV))  
CF_2016_12EPV <- as.data.frame(CF_2016_12EPV[CF_2016_12EPV!=0,])
CF_2016_12EPV <- CF_2016_12EPV[order(CF_2016_12EPV[1,])]
names(CF_2016_12EPV) <- "coef"
print(CF_2016_12EPV)

### MODELLING - GLMSPARSENET - SEPARATE2GROUPS
#USING COEFS FROM ELASTIC NET PERFORMED ON ALL 2016 INITIAL DATASET
CF_2016_12EPV <- as.matrix(coef(glasso_lasso_cox_model_643v_2016_12EPV, 
                                s = lambda_value_2016_12EPV))
sep2groups2016_12EPV <- glmSparseNet::separate2GroupsCox(as.vector(CF_2016_12EPV), 
                                                         Panglioma_2016_ALL[, rownames(CF_2016_12EPV)], 
                                                         SURV, 
                                                         plot.title = 'GLASSO+LASSO (2016_12EPV = 19 features)', legend.outside = FALSE)

sep2groups2016_12EPV

#DIVIDE INTO LOW AND HIGH GROUPS
d_12EPV_2016 <- data.frame(sep2groups2016_12EPV[["km"]][["custom.data"]])

df_low_12EPV_2016 <- d_12EPV_2016[d_12EPV_2016$group == "Low risk", ]
df_high_12EPV_2016 <- d_12EPV_2016[d_12EPV_2016$group == "High risk", ]

write.csv(df_low_12EPV_2016,"~/paraRoberta/single-omics-ds/list_LOWgroup_12EPV_2016.csv", row.names = TRUE)
write.csv(df_high_12EPV_2016,"~/paraRoberta/single-omics-ds/list_HIGHgroup_12EPV_2016.csv", row.names = TRUE)

# STATISTICS / PERCENTAGES
#AUX
low_group_patients_12EPV_2016 <- rownames(df_low_12EPV_2016)
high_group_patients_12EPV_2016 <- rownames(df_high_12EPV_2016)

#table(panglioma_clinic_single_omics_2016$histological_type)

lgg_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type %in% c('astrocytoma','oligoastrocytoma','oligodendroglioma'),]
lgg_patients <- lgg_patients$Patient_ID

gbm_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type == "untreated primary (de novo) gbm",]
gbm_patients <- gbm_patients$Patient_ID

length(lgg_patients) #494 LGG
length(gbm_patients) #149 GBM

percentagem_LGG_LOWrisk_group <- (sum(low_group_patients_12EPV_2016 %in% lgg_patients) * 100 ) / length(low_group_patients_12EPV_2016)
round(percentagem_LGG_LOWrisk_group,1)

percentagem_GBM_HIGHrisk_group <- (sum(high_group_patients_12EPV_2016 %in% gbm_patients) * 100 ) / length(high_group_patients_12EPV_2016)
round(percentagem_GBM_HIGHrisk_group,1)


sum(low_group_patients_10EPV_2016 %in% low_group_patients_12EPV_2016)
sum(high_group_patients_12EPV_2016 %in% high_group_patients_10EPV_2016)






# 15 EPV ------------------------------------------------------------------
# = 233 / 15 = 15,53
glasso_lasso_cox_model_643v_2016_15EPV <- glmnet(Panglioma_2016_ALL, 
                                                 Surv(SURV$time, SURV$status),
                                                 alpha = 1, 
                                                 family = "cox",
                                                 pmax = 19)

print(glasso_lasso_cox_model_643v_2016_15EPV)


ix_2016_15EPV           <- length(glasso_lasso_cox_model_643v_2016_15EPV$lambda)
lambda_value_2016_15EPV <- glasso_lasso_cox_model_643v_2016_15EPV$lambda[ix_2016_15EPV]
model_2016_15EPV        <- glmnet(Panglioma_2016_ALL, Surv(SURV$time, SURV$status),
                                  alpha = 1,
                                  lambda = lambda_value_2016_15EPV,
                                  family = "cox")

# extract coefficients at a single value of lambda
CF_2016_15EPV <- as.matrix(coef(glasso_lasso_cox_model_643v_2016_15EPV, s = lambda_value_2016_15EPV))  
CF_2016_15EPV <- as.data.frame(CF_2016_15EPV[CF_2016_15EPV!=0,])
CF_2016_15EPV <- CF_2016_15EPV[order(CF_2016_15EPV[1,])]
names(CF_2016_15EPV) <- "coef"
print(CF_2016_15EPV)

### MODELLING - GLMSPARSENET - SEPARATE2GROUPS
#USING COEFS FROM ELASTIC NET PERFORMED ON ALL 2016 INITIAL DATASET
CF_2016_15EPV <- as.matrix(coef(glasso_lasso_cox_model_643v_2016_15EPV, 
                                s = lambda_value_2016_15EPV))
sep2groups2016_15EPV <- glmSparseNet::separate2GroupsCox(as.vector(CF_2016_15EPV), 
                                                         Panglioma_2016_ALL[, rownames(CF_2016_15EPV)], 
                                                         SURV, 
                                                         plot.title = 'GLASSO+LASSO (2016_15EPV = 15 features)', legend.outside = FALSE)

sep2groups2016_15EPV

#DIVIDE INTO LOW AND HIGH GROUPS
d_15EPV_2016 <- data.frame(sep2groups2016_15EPV[["km"]][["custom.data"]])

df_low_15EPV_2016 <- d_15EPV_2016[d_15EPV_2016$group == "Low risk", ]
df_high_15EPV_2016 <- d_15EPV_2016[d_15EPV_2016$group == "High risk", ]

write.csv(df_low_15EPV_2016,"~/paraRoberta/single-omics-ds/list_LOWgroup_15EPV_2016.csv", row.names = TRUE)
write.csv(df_high_15EPV_2016,"~/paraRoberta/single-omics-ds/list_HIGHgroup_15EPV_2016.csv", row.names = TRUE)

# STATISTICS / PERCENTAGES
#AUX
low_group_patients_15EPV_2016 <- rownames(df_low_15EPV_2016)
high_group_patients_15EPV_2016 <- rownames(df_high_15EPV_2016)

#table(panglioma_clinic_single_omics_2016$histological_type)

lgg_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type %in% c('astrocytoma','oligoastrocytoma','oligodendroglioma'),]
lgg_patients <- lgg_patients$Patient_ID

gbm_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type == "untreated primary (de novo) gbm",]
gbm_patients <- gbm_patients$Patient_ID

length(lgg_patients) #494 LGG
length(gbm_patients) #149 GBM

percentagem_LGG_LOWrisk_group <- (sum(low_group_patients_15EPV_2016 %in% lgg_patients) * 100 ) / length(low_group_patients_15EPV_2016)
round(percentagem_LGG_LOWrisk_group,1)

percentagem_GBM_HIGHrisk_group <- (sum(high_group_patients_15EPV_2016 %in% gbm_patients) * 100 ) / length(high_group_patients_15EPV_2016)
round(percentagem_GBM_HIGHrisk_group,1)


sum(low_group_patients_10EPV_2016 %in% low_group_patients_15EPV_2016)
sum(high_group_patients_15EPV_2016 %in% high_group_patients_10EPV_2016)

sum(low_group_patients_12EPV_2016 %in% low_group_patients_15EPV_2016)
sum(high_group_patients_12EPV_2016 %in% high_group_patients_15EPV_2016)



# 10% dataset size --------------------------------------------------------
# = 643 x 0.10 = 64,3
glasso_lasso_cox_model_643v_2016_10percent <- glmnet(Panglioma_2016_ALL, 
                                                 Surv(SURV$time, SURV$status),
                                                 alpha = 1, 
                                                 family = "cox",
                                                 pmax = 95)

print(glasso_lasso_cox_model_643v_2016_10percent)


ix_2016_10percent           <- length(glasso_lasso_cox_model_643v_2016_10percent$lambda)
lambda_value_2016_10percent <- glasso_lasso_cox_model_643v_2016_10percent$lambda[ix_2016_10percent]
model_2016_10percent        <- glmnet(Panglioma_2016_ALL, Surv(SURV$time, SURV$status),
                                  alpha = 1,
                                  lambda = lambda_value_2016_10percent,
                                  family = "cox")

# extract coefficients at a single value of lambda
CF_2016_10percent <- as.matrix(coef(glasso_lasso_cox_model_643v_2016_10percent, s = lambda_value_2016_10percent))  
CF_2016_10percent <- as.data.frame(CF_2016_10percent[CF_2016_10percent!=0,])
CF_2016_10percent <- CF_2016_10percent[order(CF_2016_10percent[1,])]
names(CF_2016_10percent) <- "coef"
print(CF_2016_10percent)

### MODELLING - GLMSPARSENET - SEPARATE2GROUPS
#USING COEFS FROM ELASTIC NET PERFORMED ON ALL 2016 INITIAL DATASET
CF_2016_10percent <- as.matrix(coef(glasso_lasso_cox_model_643v_2016_10percent, 
                                s = lambda_value_2016_10percent))
sep2groups2016_10percent <- glmSparseNet::separate2GroupsCox(as.vector(CF_2016_10percent), 
                                                         Panglioma_2016_ALL[, rownames(CF_2016_10percent)], 
                                                         SURV, 
                                                         plot.title = 'GLASSO+LASSO (2016_10percent = 64 features)', legend.outside = FALSE)

sep2groups2016_10percent

#DIVIDE INTO LOW AND HIGH GROUPS
d_10percent_2016 <- data.frame(sep2groups2016_10percent[["km"]][["custom.data"]])

df_low_10percent_2016 <- d_10percent_2016[d_10percent_2016$group == "Low risk", ]
df_high_10percent_2016 <- d_10percent_2016[d_10percent_2016$group == "High risk", ]

write.csv(df_low_10percent_2016,"~/paraRoberta/single-omics-ds/list_LOWgroup_10percent_2016.csv", row.names = TRUE)
write.csv(df_high_10percent_2016,"~/paraRoberta/single-omics-ds/list_HIGHgroup_10percent_2016.csv", row.names = TRUE)

# STATISTICS / PERCENTAGES
#AUX
low_group_patients_10percent_2016 <- rownames(df_low_10percent_2016)
high_group_patients_10percent_2016 <- rownames(df_high_10percent_2016)

#table(panglioma_clinic_single_omics_2016$histological_type)

lgg_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type %in% c('astrocytoma','oligoastrocytoma','oligodendroglioma'),]
lgg_patients <- lgg_patients$Patient_ID

gbm_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type == "untreated primary (de novo) gbm",]
gbm_patients <- gbm_patients$Patient_ID

length(lgg_patients) #494 LGG
length(gbm_patients) #149 GBM

percentagem_LGG_LOWrisk_group <- (sum(low_group_patients_10percent_2016 %in% lgg_patients) * 100 ) / length(low_group_patients_10percent_2016)
round(percentagem_LGG_LOWrisk_group,1)

percentagem_GBM_HIGHrisk_group <- (sum(high_group_patients_10percent_2016 %in% gbm_patients) * 100 ) / length(high_group_patients_10percent_2016)
round(percentagem_GBM_HIGHrisk_group,1)


sum(low_group_patients_10EPV_2016 %in% low_group_patients_10percent_2016)
sum(high_group_patients_10EPV_2016 %in% high_group_patients_10percent_2016)

sum(low_group_patients_12EPV_2016 %in% low_group_patients_10percent_2016)
sum(high_group_patients_12EPV_2016 %in% high_group_patients_10percent_2016)

sum(low_group_patients_15EPV_2016 %in% low_group_patients_10percent_2016)
sum(high_group_patients_15EPV_2016 %in% high_group_patients_10percent_2016)




# 15% dataset size --------------------------------------------------------
# = 643 x 0.15 = 96,45

glasso_lasso_cox_model_643v_2016_15percent <- glmnet(Panglioma_2016_ALL, 
                                                     Surv(SURV$time, SURV$status),
                                                     alpha = 1, 
                                                     family = "cox",
                                                     pmax = 130)

print(glasso_lasso_cox_model_643v_2016_15percent)


ix_2016_15percent           <- length(glasso_lasso_cox_model_643v_2016_15percent$lambda)
lambda_value_2016_15percent <- glasso_lasso_cox_model_643v_2016_15percent$lambda[ix_2016_15percent]
model_2016_15percent        <- glmnet(Panglioma_2016_ALL, Surv(SURV$time, SURV$status),
                                      alpha = 1,
                                      lambda = lambda_value_2016_15percent,
                                      family = "cox")

# extract coefficients at a single value of lambda
CF_2016_15percent <- as.matrix(coef(glasso_lasso_cox_model_643v_2016_15percent, s = lambda_value_2016_15percent))  
CF_2016_15percent <- as.data.frame(CF_2016_15percent[CF_2016_15percent!=0,])
CF_2016_15percent <- CF_2016_15percent[order(CF_2016_15percent[1,])]
names(CF_2016_15percent) <- "coef"
print(CF_2016_15percent)

### MODELLING - GLMSPARSENET - SEPARATE2GROUPS
#USING COEFS FROM ELASTIC NET PERFORMED ON ALL 2016 INITIAL DATASET
CF_2016_15percent <- as.matrix(coef(glasso_lasso_cox_model_643v_2016_15percent, 
                                    s = lambda_value_2016_15percent))
sep2groups2016_15percent <- glmSparseNet::separate2GroupsCox(as.vector(CF_2016_15percent), 
                                                             Panglioma_2016_ALL[, rownames(CF_2016_15percent)], 
                                                             SURV, 
                                                             plot.title = 'GLASSO+LASSO (2016_15percent = 97 features)', legend.outside = FALSE)

sep2groups2016_15percent

#DIVIDE INTO LOW AND HIGH GROUPS
d_15percent_2016 <- data.frame(sep2groups2016_15percent[["km"]][["custom.data"]])

df_low_15percent_2016 <- d_15percent_2016[d_15percent_2016$group == "Low risk", ]
df_high_15percent_2016 <- d_15percent_2016[d_15percent_2016$group == "High risk", ]

write.csv(df_low_15percent_2016,"~/paraRoberta/single-omics-ds/list_LOWgroup_15percent_2016.csv", row.names = TRUE)
write.csv(df_high_15percent_2016,"~/paraRoberta/single-omics-ds/list_HIGHgroup_15percent_2016.csv", row.names = TRUE)

# STATISTICS / PERCENTAGES
#AUX
low_group_patients_15percent_2016 <- rownames(df_low_15percent_2016)
high_group_patients_15percent_2016 <- rownames(df_high_15percent_2016)

#table(panglioma_clinic_single_omics_2016$histological_type)

lgg_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type %in% c('astrocytoma','oligoastrocytoma','oligodendroglioma'),]
lgg_patients <- lgg_patients$Patient_ID

gbm_patients <- panglioma_clinic_single_omics_2016[panglioma_clinic_single_omics_2016$histological_type == "untreated primary (de novo) gbm",]
gbm_patients <- gbm_patients$Patient_ID

length(lgg_patients) #494 LGG
length(gbm_patients) #149 GBM

percentagem_LGG_LOWrisk_group <- (sum(low_group_patients_15percent_2016 %in% lgg_patients) * 100 ) / length(low_group_patients_15percent_2016)
round(percentagem_LGG_LOWrisk_group,1)

percentagem_GBM_HIGHrisk_group <- (sum(high_group_patients_15percent_2016 %in% gbm_patients) * 100 ) / length(high_group_patients_15percent_2016)
round(percentagem_GBM_HIGHrisk_group,1)


sum(low_group_patients_10EPV_2016 %in% low_group_patients_15percent_2016)
sum(high_group_patients_10EPV_2016 %in% high_group_patients_15percent_2016)

sum(low_group_patients_12EPV_2016 %in% low_group_patients_15percent_2016)
sum(high_group_patients_12EPV_2016 %in% high_group_patients_15percent_2016)

sum(low_group_patients_15EPV_2016 %in% low_group_patients_15percent_2016)
sum(high_group_patients_15EPV_2016 %in% high_group_patients_15percent_2016)



























# OLD ---------------------------------------------------------------------



# OLD
### MODELLING
# 0. Model preparation

# SELECT BALANCED FOLDS FOR CROSS-VALIDATION
foldid <- balanced.cv.folds(!!SURV$status, nfolds = 10)$output #("!!" converts 1 and 0 to TRUE and FALSE)

# List that will store all selected genes
selected.genes <- list()

##########################################################################
#0. Stratified Regularized Cox Regression Model
#https://glmnet.stanford.edu/articles/Coxnet.html
strata <- rep(1:5, length.out = 643)
y <- Surv(SURV$time, SURV$status)
y2 <- stratifySurv(y, strata)
str(y2[1:6])
#stratifySurv returns an object of class stratifySurv. We can then pass this
#stratifySurv object as the response to a glmnet call. glmnet will fit a
#stratified Cox model if it detects that the response has class stratifySurv
x   <- as.matrix(Pan_glioma_i_2016)
######fit <- glmnet(x, y2, family = "cox")
plot(fit)
######save(fit, file = "~/paraRoberta/fit.RData")
load(file = "~/paraRoberta/fit.RData")
#This stratifySurv object can also be passed to cv.glmnet to fit stratified Cox models with cross-validation:
######cv.fit <- cv.glmnet(x, y2, family = "cox", nfolds = 5)
######save(cv.fit, file = "~/paraRoberta/cv.fit.RData")
load(file = "~/paraRoberta/cv.fit.RData")

plot(cv.fit)

plot(survival::survfit(cv.fit, s = cv.fit$lambda.min, x = x, y = y2))
#Plotting survival curves
#Fitting a regularized Cox model using glmnet with family = "cox" returns an object of class "coxnet".
#Class "coxnet" objects have a survfit method which allows the user to visualize the survival curves from the model.
#In addition to the "coxnet" object, the user must pass the x and y objects used to fit the model
#(for computation of the baseline hazard), as well as the lambda value for which the survival curve is desired:


#Coefficients of selected model from Cross-Validation

#Taking the best model described by lambda.min

coefs.v <- coef(cv.fit, s = 'lambda.min')[,1] %>% { .[. != 0]}
coefs.v %>% { 
  data.frame(gene.name   = names(.),
             coefficient = .,
             stringsAsFactors = FALSE)
} %>%
  arrange(gene.name) %>%
  knitr::kable()


### MODELLING - GLMSPARSENET - SEPARATE2GROUPS

#USING COEFS FROM ELASTIC NET PERFORMED ON ALL 2016 INITIAL DATASET
sep2groups2016 <- glmSparseNet::separate2GroupsCox(as.vector(coefs.v), 
                                 Pan_glioma_i_2016[, names(coefs.v)], 
                                 SURV, 
                   plot.title = 'GLASSO all 2016 variables', legend.outside = FALSE)

sep2groups2016

#DIVIDE INTO LOW AND HIGH GROUPS
d <- data.frame(sep2groups2016[["km"]][["custom.data"]])

df_low <- d[d$group == "Low risk", ]
df_high <- d[d$group == "High risk", ]



#IDENTIFYING PATIENTS OF EACH GROUP
d <- data.frame(time = sep2groups2016$km$time,
                n.risk = sep2groups2016$km$n.risk,
                n.event = sep2groups2016$km$n.event,
                n.censor = sep2groups2016$km$n.censor,
                surv = sep2groups2016$km$surv,
                upper = sep2groups2016$km$upper,
                lower = sep2groups2016$km$lower
)
head(d)

aux <- sep2groups2016$plot$plot$data
sep2groups2016$plot$plot$plot_env$df
sep2groups2016$plot$plot$plot_env$fit
sep2groups2016$plot$data.survplot$strata


#USING COEFS FROM ELASTIC NET PERFORMED ON ALL 2016 INITIAL DATASET




## 1. Elastic Net (classic penalization method) ----------------------------------------------------------
# - Uses regular glmnet model as simple baseline
# Does k-fold cross-validation for glmnet, produces a plot, and returns a value for lambda
Pan_glioma_i_2016 <- as.matrix(Pan_glioma_i_2016)
cv.glm <- cv.glmnet(Pan_glioma_i_2016, 
                    Surv(SURV$time, SURV$status), 
                    family = 'cox',
                    foldid = foldid)

#save object
saveRDS(cv.glm, file = "cv_elasticnet.rds")
# Restore the object
#cv.glm <- readRDS(file = "cv_elasticnet.rds")

glmSparseNet::separate2GroupsCox(as.vector(coef(cv.glm, s = 'lambda.min')[,1]), 
                                 as.matrix(Pan_glioma_i_2016), SURV, 
                                 plot.title = 'Full dataset', 
                                 legend.outside = FALSE)

selected.genes[['GLMnet']] <- coef(cv.glm, s = 'lambda.min')[,1] %>% 
  { .[. != 0] } %>%
  names %>% 
  geneNames %>% 
  { .[['external_gene_name']]} 

selected.genes[['GLMnet']]




#comparing nested models using Likelihood Ratio Test (LRT) (can we drop AACS? without making the model significantly worse)
#Estimating survival as a function of "AACS" and "ABCA10" measures.
cox.mod1 <- coxph(Surv(time, status) ~ AACS + ABCA10, data = Pan_glioma_i_2016_and_SURV)
cox.mod2 <- coxph(Surv(time, status) ~ ABCA10,        data = Pan_glioma_i_2016_and_SURV)
# LRT to compare the two models
anova(cox.mod2, cox.mod1, test="LRT")
# P-value is large, thus, there is not a statistically significant difference between the two models.
# We can drop the "AACS" without loss of predictive power

summary(cox.mod1)
summary(cox.mod2)

#Create survival curves
cox_fit1 <- survfit(cox.mod1)
cox_fit2 <- survfit(cox.mod2)
autoplot(cox_fit1, xlab = "Time (days)", ylab = "Survival (%)", main = "1 - Cox Proportional Hazards Model Plot: TCGA-PanGlioma 2016 cases (n=643)")
autoplot(cox_fit2, xlab = "Time (days)", ylab = "Survival (%)", main = "1 - Cox Proportional Hazards Model Plot: TCGA-PanGlioma 2016 cases (n=643)")



##########################################################
### CHECKING COX PROPORTIONAL HAZARD MODEL ASSUMPTIONS ###
##########################################################


# 1 - Checking the linearity assumptions:
#We assume that the relationship between any of the numeric x variables and the log-hazard is linear
#How? As with any other model. We look into the residual plot
# x axis: plot the predicted values, y axis: the residuals
#(for survival analysis: the Martingale residuals, but could also be the deviance residuals)



plot(predict(cox.mod1), residuals(cox.mod1, type = "martingale"),
     xlab = "fitted values", ylab = "Martingale residuals", 
     main = "Residual Plot", las = 1) #las = rotates the values on the y axis
#add a line ax y=residual=0
abline(h=0) #horizontal line
#fit a smoother through the points
lines(smooth.spline(predict(cox.mod1), residuals(cox.mod1, type = "martingale")), col="red")

#DRAW KAPLAN MEYER SURVIVAL CURVES
#Calculate log rank test
library(glmSparseNet)
separate2GroupsCox(chosen.btas = , Pan_glioma_f_2016, SURV_reordered, plot.title = "all 16K")

#OLD
# DATA LOADING
#Load RData
#2016 SINGLE OMICS
#load("~/paraRoberta/single-omics-ds/single-omics_starting_ds_2016.RData") #Pan_glioma_i & Pan_glioma_norm_i #643 patients
#load("~/paraRoberta/single-omics-ds/single-omics_final_ds_2016.RData")    #Pan_glioma_f & Pan_glioma_norm_f #643 patients 
#Pan_glioma_i_2016 <- Pan_glioma_i[-(dim(Pan_glioma_i)[2])] #remove last column = names
#Pan_glioma_f_2016 <- Pan_glioma_f
#rm(Pan_glioma_i, Pan_glioma_f)


#2021 SINGLE OMICS
#load("~/paraRoberta/single-omics-ds/single-omics_starting_ds_2021.RData") #Pan_glioma_i & Pan_glioma_norm_i #583 patients
#load("~/paraRoberta/single-omics-ds/single-omics_final_ds_2021.RData")    #Pan_glioma_f & Pan_glioma_norm_f #583 patients
#Pan_glioma_i_2021 <- Pan_glioma_i[-(dim(Pan_glioma_i)[2])] #remove last column = names
#Pan_glioma_f_2021 <- Pan_glioma_f
#rm(Pan_glioma_i, Pan_glioma_f, Pan_glioma_norm_i, Pan_glioma_norm_f)




###############################################################################
####################### UNIVARIATE COX REGRESSION #############################
###############################################################################
#To apply the univariate coxph function to multiple covariates at once, type this:
#covariates_2016_glasso_selected
covariates <- colnames(Panglioma_2016_ALL)
covariates2 <- colnames(Panglioma_2016_ALL)

#correcting "-" character from genes' names bc it gives problems for the formula of the model:
colnames(Pan_glioma_f_2016_and_SURV) <- sub("-", "", colnames(Pan_glioma_f_2016_and_SURV), fixed = TRUE)
covariates2 <- sub("-", "", covariates2, fixed = TRUE)

#covariates <- c("A2BP1", "ADAM11")
#covariates <- c("age", "sex",  "ph.karno", "ph.ecog", "wt.loss")
univ_formulas <- sapply(covariates2,
                        function(x) as.formula(paste('Surv(time, status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = Pan_glioma_f_2016_and_SURV)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(apply(res, 2, unlist)) #!!!!!!!
res$p.value <- as.numeric(res$p.value)

#RESULTS / FILTERING
res_significant005   <- res[res$p.value <= 0.05, ]
res_significant001   <- res[res$p.value <= 0.01, ]
res_significant0001  <- res[res$p.value <= 0.001, ]
res_significant00001 <- res[res$p.value <= 0.0001, ]
res_significant000001 <- res[res$p.value <= 0.00001, ]
res_significant0000000001 <- res[res$p.value <= 0.000000001, ]


#The output above shows the regression beta coefficients, the effect sizes 
#(given as hazard ratios) and statistical significance for each of the variables 
#in relation to overall survival. Each factor is assessed through separate univariate Cox regressions.



###############################################################################
####################### MULTIVARIATE COX REGRESSION ANALYSIS ##################
###############################################################################
#Describe how the factors jointly impact on survival.
#we perform a multivariate Cox regression analysis.
#We remove non-statistically significant variables (from Univariate Cox Regression Analysis' results)
#in the following multivariate analysis.
#We include the X factors into the multivariate model.

STRING <- paste0(rownames(res_significant0000000001), collapse = " + ")

res.cox0000000001 <- coxph(Surv(time, status) ~ ABI3 + ACTL6B + ACVRL1 + ADAP2 + ADORA3 + AGAP4 + ALOX5 + ANG + ANXA2P1 + ANXA2P2 + ANXA2 + ARHGAP11A + ARHGAP9 + ARHGDIB + ARL11 + ARPC1B + ASF1B + ASPM + ATP1A3 + ATP5J2 + AURKA + AURKB + BIRC5 + BIRC6 + BRIP1 + BTK + BUB1 + C16orf54 + C17orf60 + C19orf56 + C1QA + C1QB + C1QC + C1R + C1S + C1orf151 + C1orf38 + C20orf177 + C3 + C7orf59 + C9orf102 + CCL5 + CCNA2 + CCNB2 + CCR1 + CCR5 + CD14 + CD2 + CD300A + CD300C + CD33 + CD37 + CD3E + CD4 + CD53 + CD68 + CD74 + CD86 + CDC20 + CDC25C + CDC45 + CDC6 + CDCA2 + CDCA5 + CDCA8 + CDK1 + CELF3 + CELF5 + CENPA + CENPE + CENPF + CENPK + CEP55 + CKAP2L + CKMT1B + CLEC18A + CLEC18B + CLEC7A + CLIC1 + CLSPN + CMPK2 + CMTM7 + COL1A1 + COL1A2 + COL4A1 + COL4A2 + COX6B1 + CSDAP1 + CSDA + CSF2RB + CSF3R + CSRNP3 + CTSS + CXCR6 + CYBA + CYTH4 + DAPP1 + DEPDC1B + DEPDC1 + DHRS4 + DLGAP5 + DTL + DTX3L + E2F2 + ECT2 + EEF1A1P9 + EEF1G + EIF3L + EIF5AL1 + EME1 + EPR1 + EPSTI1 + ERCC6L + ESAM + ESCO2 + ESPL1 + EVI2B + EXO1 + EXOC6B + FAM111B + FAM115C + FAM123C + FAM64A + FAM72B + FAM72D + FANCA + FANCD2 + FANCI + FCER1G + FCGR1A + FCGR2A + FCGR3A + FERMT3 + FOXM1 + FOXO3B + FOXO3 + FPR1 + GMFG + GNA15 + GNGT2 + GPR65 + GPSM3 + GRAP + GSG2 + GTF2A1 + GTF2IRD2 + GTSE1 + HAVCR2 + HCK + HCLS1 + HJURP + HLAA + HLADMA + HLADMB + HLADPA1 + HLADPB1 + HLADQA1 + HLADQB1 + HLADRA + HLADRB1 + HLAH + HMMR + HMP19 + HNRNPA1L2 + HNRNPA1 + IFI30 + IFI6 + IL10RA + IL18 + INA + IQGAP3 + ITGB2 + KIAA0101 + KIF11 + KIF14 + KIF20A + KIF23 + KIF2C + KIF4A + KIF4B + KIFC1 + KLHL26 + KPNA5 + LAIR1 + LAPTM5 + LAT2 + LCK + LCP1 + LGALS9 + LOC100233209 + LOC154761 + LOC541471 + LOC606724 + LOC653566 + LRRC25 + MAPK8IP2 + MELK + MLF1IP + MRPS15 + MS4A4A + MS4A6A + MSR1 + MX1 + MYBL2 + MYH6 + MYH7 + MYL12A + MYO1F + MYO9A + MYST3 + NACA + NCAPG + NCAPH + NCF4 + NCOA2 + NCRNA00152 + NDC80 + NEIL3 + NEK2 + NSUN5P1 + NUDT16P1 + NUF2 + NUSAP1 + OAS2 + OLIG1 + OLIG2 + ORC1L + ORC6L + OSCAR + PABPC3 + PARP9 + PARVG + PBK + PDIA3P + PDIA3 + PHF5A + PIK3AP1 + PIK3CG + PKMYT1 + PLEKHM3 + PLEK + PLK1 + PLK4 + PNMAL2 + POC1A + POLQ + POLR2I + POM121C + PRC1 + PSD + PSMB8 + PSMB9 + PTGS1 + PTPRC + PYCARD + RAD51 + RAD54L + RBM47 + RECQL4 + RGS19 + RIPK3 + RNASE2 + RNASE4 + RNASE6 + RNASET2 + RPL10 + RPL12 + RPL13A + RPL18A + RPL3 + RPL4 + RPL6 + RPL7A + RPL7 + RPS12 + RPS24 + RPS25 + RPS27A + RPS6KA1 + RRM2 + RUNDC3A + S100A11 + S100A8 + S100A9 + SAMD8 + SAMSN1 + SASH3 + SCRT1 + SGOL1 + SGOL2 + SHCBP1 + SIGLEC7 + SIGLEC9 + SIPA1 + SKA1 + SKA3 + SLA + SLC15A3 + SLC2A3 + SLC7A7 + SMG1 + SNAP91 + SOCS7 + SPAG5 + SPC24 + SPC25 + SPCS2 + SPI1 + SRGN + STX1B + STXBP2 + SUMO1P3 + SVOP + SYNGR2 + TAOK1 + TBCB + TFEC + TK1 + TLR1 + TLR2 + TMEM106A + TMEM176A + TMEM176B + TOP2A + TPX2 + TREM2 + TROAP + TTBK2 + TTC30B + TTK + TYROBP + UBE2C + UBE2MP1 + UBL5 + UNC13A + UQCR11 + UQCRHL + UQCRQ + USP34 + VAMP8 + VSIG4 + WDR62 + ZKSCAN1 + ZNF778 + ZNHIT1 + ZWINT
                      , iter.max = 500, data =  Pan_glioma_f_2016_and_SURV)
summary(res.cox0000000001)

#coxph.fit(Panglioma_2016_ALL, SURV, strata=NULL)



# Plot the baseline survival function
ggsurvplot(survfit(res.cox0000000001), color = "#2E9FDF", ggtheme = theme_minimal())
