# # # # # # # # # # # # # # # # # # # #
#                                    #
#             FIGURES                 ----
#     Space-for-time manuscript      #
#                                    #
# # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # #
# CONTENT                  ----
# AIM : contains code to create all figures
# The data paths need to be adapted to the directory situation. If e.g. the recommendations
# of folder Structure in 1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/README.txt were followed:
# "1.Dataset creation/step 1 - Prepare_input_data_plotwise_pairwise/data/OutputData/"


# # # # # # # # # # # # # #
# 0 - REQUIREMENTS              ----
#
# packages
library(ggplot2)
library(patchwork)
library(vegan)

# # # # #
# 0.a. - DATA ----
#data
#set working directory to folder "4. Plotting results"
#pairwise datasets
load("./data/InputData/pwise_time_plants.RData")
load("./data/InputData/pwise_time_herb.RData")
load("./data/InputData/pwise_time_pred.RData")

load("./data/InputData/pwise_space_plants.RData")
load("./data/InputData/pwise_space_herb.RData")
load("./data/InputData/pwise_space_pred.RData")


#GDM splines
#alpha
load("./data/InputData/splines_alpha_plant_sc_space.RData")
load("./data/InputData/splines_alpha_plant_sc_time.RData")

load("./data/InputData/splines_alpha_pred_sc_space.RData")
load("./data/InputData/splines_alpha_pred_sc_time.RData")

load("./data/InputData/splines_alpha_herb_sc_space.RData")
load("./data/InputData/splines_alpha_herb_sc_time.RData")

#beta
load("./data/InputData/splines_beta_plant_sc_space.RData")
load("./data/InputData/splines_beta_plant_sc_time.RData")

load("./data/InputData/splines_beta_pred_sc_space.RData")
load("./data/InputData/splines_beta_pred_sc_time.RData")

load("./data/InputData/splines_beta_herb_sc_space.RData")
load("./data/InputData/splines_beta_herb_sc_time.RData")

# 0.b  -  produce merged datasets for later plots ####

# # # # code how to produce div_temp_fig, and div_space_fig.RData # # # # # # # # # # # # # # #  # # # #
# dataset having alpha and beta diversity of all organisms together
# # # # # # # # # # # # # # #  # # # #

#in case it does not work, it might be that the vector of standardised differences
#in alpha diversity - da0/1/2/3/4st got lost

#so check if $da0/1/2/3/4st, $dha0/1/2/3/4st, $dpa0/1/2/3/4st is there, 
# if missing run code in folder 1. dataset creation/1.Prepare input data...

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#temporal diva/group
Bsim_time<-c(pwise_time_plants$bsim, #plants
             pwise_time_herb$hbsim,#herbivores
             pwise_time_pred$pbsim)#sec. consumers

B0_time<-c(pwise_time_plants$cqn0dis, 
           pwise_time_herb$hcqn0dis,
           pwise_time_pred$pcqn0dis) 

B1_time<-c(pwise_time_plants$cqn1dis,
           pwise_time_herb$hcqn1dis,
           pwise_time_pred$pcqn1dis) 

B2_time<-c(pwise_time_plants$cqn2dis,
           pwise_time_herb$hcqn2dis,
           pwise_time_pred$pcqn2dis)

B3_time<-c(pwise_time_plants$cqn3dis,
           pwise_time_herb$hcqn3dis,
           pwise_time_pred$pcqn3dis) 

B4_time<-c(pwise_time_plants$cqn4dis,
           pwise_time_herb$hcqn4dis,
           pwise_time_pred$pcqn4dis) 

A0_time<-c(pwise_time_plants$da0st, 
           pwise_time_herb$dha0st,
           pwise_time_pred$dpa0st) 

A1_time<-c(pwise_time_plants$da1st, 
           pwise_time_herb$dha1st,
           pwise_time_pred$dpa1st) 

A2_time<-c(pwise_time_plants$da2st, 
           pwise_time_herb$dha2st,
           pwise_time_pred$dpa2st) 

A3_time<-c(pwise_time_plants$da3st, 
           pwise_time_herb$dha3st,
           pwise_time_pred$dpa3st) 

A4_time<-c(pwise_time_plants$da4st, 
           pwise_time_herb$dha4st,
           pwise_time_pred$dpa4st) 

EP_time<-c(pwise_time_plants$EP, pwise_time_herb$EP, 
           pwise_time_pred$EP) 

dYR_time<-c(pwise_time_plants$dYR, pwise_time_herb$dYR, 
            pwise_time_pred$dYR) 

mLUIraw<- c(pwise_time_plants$mLUI_raw, pwise_time_herb$mLUI_raw, 
            pwise_time_pred$mLUI_raw)

mLUIscaled<- c(pwise_time_plants$mLUI, pwise_time_herb$mLUI, 
               pwise_time_pred$mLUI) 

dLUI<- c(pwise_time_plants$dLUI, pwise_time_herb$dLUI, 
         pwise_time_pred$dLUI) 

mMOWraw<- c(pwise_time_plants$mMOW_raw, pwise_time_herb$mMOW_raw, 
            pwise_time_pred$mMOW_raw) 

mMOWscaled<- c(pwise_time_plants$mMOW, pwise_time_herb$mMOW, 
               pwise_time_pred$mMOW) 

dMOW<- c(pwise_time_plants$dMOW, pwise_time_herb$dMOW, 
         pwise_time_pred$dMOW) 

mGRAraw<- c(pwise_time_plants$mGRA_raw, pwise_time_herb$mGRA_raw, 
            pwise_time_pred$mGRA_raw)

mGRAscaled<- c(pwise_time_plants$mGRA, pwise_time_herb$mGRA, 
               pwise_time_pred$mGRA) 

dGRA<- c(pwise_time_plants$dGRA, pwise_time_herb$dGRA, 
         pwise_time_pred$dGRA)

mFERraw<- c(pwise_time_plantss$mFER_raw, pwise_time_herb$mFER_raw, 
            pwise_time_pred$mFER_raw)

mFERscaled<- c(pwise_time_plants$mFER, pwise_time_herb$mFER, 
               pwise_time_pred$mFER) 

dFER<- c(pwise_time_plants$dFER, pwise_time_herb$dFER, 
         pwise_time_pred$dFER) 

group<- rep(c("plant", "herbivore", "secondary consumer"),
            times= c(nrow(pwise_time_plants),
                     nrow(pwise_time_herb),nrow(pwise_time_pred)))

div_temp<- data.frame(Bsim_time, B0_time, B1_time
                      ,B2_time, B3_time, B4_time,
                      A0_time, A1_time
                      ,A2_time, A3_time, A4_time,
                      EP_time, dYR_time, 
                      mLUIraw,mLUIscaled, dLUI,
                      mMOWraw, mMOWscaled, dMOW,
                      mGRAraw, mGRAscaled, dGRA,
                      mFERraw, mFERscaled, dFER,
                      group)

div_temp$group<- as.factor(div_temp$group)
div_temp$group<- factor(div_temp$group, levels =
                          c("plant","herbivore", "secondary consumer"))


#spatial div/group

Bsim_space<-c(pwise_space_plants$bsim,
              pwise_space_herb$hbsim, 
              pwise_space_pred$pbsim) 

B0_space<-c(pwise_space_plants$cqn0dis, 
            pwise_space_herb$hcqn0dis, 
            pwise_space_pred$pcqn0dis) 

B1_space<-c(pwise_space_plants$cqn1dis, 
            pwise_space_herb$hcqn1dis,
            pwise_space_pred$pcqn1dis) 

B2_space<-c(pwise_space_plants$cqn2dis, 
            pwise_space_herb$hcqn2dis,
            pwise_space_pred$pcqn2dis) 

B3_space<-c(pwise_space_plants$cqn3dis, 
            pwise_space_herb$hcqn3dis,
            pwise_space_pred$pcqn3dis) 

B4_space<-c(pwise_space_plants$cqn4dis, 
            pwise_space_herb$hcqn4dis,
            pwise_space_pred$pcqn4dis) 
A0_space<-c(pwise_space_plants$da0st, 
            pwise_space_herb$dha0st,
            pwise_space_pred$dpa0st) 

A1_space<-c(pwise_space_plants$da1st, 
            pwise_space_herb$dha1st,
            pwise_space_pred$dpa1st) 

A2_space<-c(pwise_space_plants$da2st, 
            pwise_space_herb$dha2st,
            pwise_space_pred$dpa2st) 

A3_space<-c(pwise_space_plants$da3st, 
            pwise_space_herb$dha3st,
            pwise_space_pred$dpa3st) 

A4_space<-c(pwise_space_plants$da4st, 
            pwise_space_herb$dha4st,
            pwise_space_pred$dpa4st) 

Dist_space<- c(pwise_space_plants$geo.dist,
               pwise_space_herb$geo.dist, 
               pwise_space_pred$geo.dist) 

YR_space<-c(pwise_space_plants$YR, 
            pwise_space_herb$YR, 
            pwise_space_pred$YR) 



mLUIraw<- c(pwise_space_plants$mLUI_raw, pwise_space_herb$mLUI_raw, 
            pwise_space_pred$mLUI_raw) 
mLUIscaled<- c(pwise_space_plants$mLUI, pwise_space_herb$mLUI, 
               pwise_space_pred$mLUI) 
dLUI<- c(pwise_space_plants$dLUI, pwise_space_herb$dLUI, 
         pwise_space_pred$dLUI) 

mMOWraw<- c(pwise_space_plants$mMOW_raw, pwise_space_herb$mMOW_raw, 
            pwise_space_pred$mMOW_raw) 
mMOWscaled<- c(pwise_space_plants$mMOW, pwise_space_herb$mMOW, 
               pwise_space_pred$mMOW) 
dMOW<- c(pwise_space_plants$dMOW, pwise_space_herb$dMOW, 
         pwise_space_pred$dMOW) 

mGRAraw<- c(pwise_space_plants$mGRA_raw, pwise_space_herb$mGRA_raw, 
            pwise_space_pred$mGRA_raw) 
mGRAscaled<- c(pwise_space_plants$mGRA, pwise_space_herb$mGRA, 
               pwise_space_pred$mGRA) 
dGRA<- c(pwise_space_plants$dGRA, pwise_space_herb$dGRA, 
         pwise_space_pred$dGRA)

mFERraw<- c(pwise_space_plants$mFER_raw, pwise_space_herb$mFER_raw, 
            pwise_space_pred$mFER_raw) 
mFERscaled<- c(pwise_space_plants$mFER, pwise_space_herb$mFER, 
               pwise_space_pred$mFER) 
dFER<- c(pwise_space_plants$dFER, pwise_space_herb$dFER, 
         pwise_space_pred$dFER) 
group<- rep(c("plant",  "herbivore", "secondary consumer"),
            times=c(nrow(pwise_space), nrow(pwise_space_in),nrow(pwise_space_in)))

div_space<- data.frame(Bsim_space, B0_space, B1_space
                       ,B2_space, B3_space, B4_space,
                       A0_space, A1_space, A2_space,
                       A3_space, A4_space,
                       YR_space, Dist_space, 
                       mLUIraw,mLUIscaled, dLUI,
                       mMOWraw, mMOWscaled, dMOW,
                       mGRAraw, mGRAscaled, dGRA,
                       mFERraw, mFERscaled, dFER,
                       group)

div_space$group<- as.factor(div_space$group)
div_space$group<- factor(div_space$group, levels =
                           c("plant", "herbivore",  "secondary consumer"))

reg<- c(paste(pwise_time_plants$REG), paste(pwise_time_herb$REG), 
        paste(pwise_time_pred$REG))
div_temp$reg<- as.factor(reg)
reg<- c(paste(pwise_space_plants$dREG), paste(pwise_space_herb$dREG), 
        paste(pwise_space_pred$dREG))
div_space$reg<- as.factor(reg)

#set working directory to 4. plotting results
save(div_temp, file="/data/InputData/div_temp_fig.RData")

save(div_space, file="/data/InputData/div_space_fig.RData")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # #
# 1 - GDM plots              ----
#
# merge datasets to make a GDM plot for all groups together

alpha_space<- rbind(splines_alpha_plant_space_sc, splines_alpha_herb_sc_space, splines_alpha_pred_sc_space)
alpha_space$organism<- rep(c("plants", "herbivores", "sec.consumers"), each=1000)

alpha_time<- rbind(splines_alpha_plant_time_sc, splines_alpha_herb_sc_time, splines_alpha_pred_sc_time)
alpha_time$organism<- rep(c("plants", "herbivores", "sec.consumers"), each=1000)

beta_space<- rbind(splines_plants_space_sc, splines_her_space_sc, splines_pred_space_sc)
beta_space$organism<- rep(c("plants", "herbivores", "sec.consumers"), each=1200)

beta_time<- rbind(splines_plants_yall_sc, splines_her_yall_sc, splines_pred_yall_sc)
beta_time$organism<- rep(c("plants", "herbivores", "sec.consumers"), each=1200)

# checking min and max for plotting
min(beta_space$LUIx) # -1.675092
max(beta_space$LUIx) # 2.357676
max(beta_space$LUIy[div_space$beta_type=="bsim"]) # 0.65

min(alpha_space$LUIx[alpha_space$alpha_type=="a0"]) # -1.675092
max(alpha_space$LUIx[alpha_space$alpha_type=="a0"]) # 2.357676
max(alpha_space$LUIy[alpha_space$alpha_type=="a0"]) # 0.4705429


# # # # # # # # # # # # # #
# 1.1 - ALPHA DIVERSITY              ----
#
# alpha richness

# # # # #
# 1.1.a. - MAIN FIGURES ----
#
# main figures
# alpha - richness

alpha_time$organism<- factor(alpha_time$organism, levels = c("plants", "herbivores", "sec.consumers"))
p2<- ggplot(alpha_time[alpha_time$alpha_type=="a0",], aes(x=LUIx, y=LUIy, colour=organism))+
  geom_line(size=1.5)+
  xlim(-2.5, 2.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))+
  xlab("LUI (scaled to EP average)")+
  ylab("Effect of LUI differences on alpha-diversity")+
  ggtitle("b - temporal land use intensity(LUI)\neffects on alpha diversity")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  
  theme_classic()+
  theme(text=element_text(size = 20))
p2

alpha_space$organism<- factor(alpha_space$organism, levels = c("plants", "herbivores", "sec.consumers"))
p1<- ggplot(alpha_space[alpha_space$alpha_type=="a0",], aes(x=LUIx, y=LUIy, colour=organism))+
  geom_line(size=1.5)+
  xlim(-2.5, 2.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))+
  xlab("LUI (scaled to year average)")+
  ylab("Effect of LUI differences on alpha-diversity")+
  ggtitle("a - spatial land use intensity (LUI)\neffects on alpha diversity")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))
p1


p1+p2+plot_layout(guides = "collect")


# # # # #
# 1.1.b. - SUPPLEMENTARY FIGURES ----
#
# Supplementary figures
# plotting as png with 1500*1500 pixels


# # # # #
# 1.1.b.a - plants space ----
#
#plants
#space
p1<- ggplot(splines_alpha_plant_space_sc, aes(x=LUIx, y=LUIy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("LUI (scaled to year average)")+
  ylab("Effect of LUI differences on alpha-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p1

p2<- ggplot(splines_alpha_plant_space_sc, aes(x=MOWx, y=MOWy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("MOW (scaled to year average)")+
  ylab("Effect of MOW differences on alpha-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_alpha_plant_space_sc, aes(x=GRAx, y=GRAy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("GRA (scaled to year average)")+
  ylab("Effect of GRA differences on alpha-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p3


p4<- ggplot(splines_alpha_plant_space_sc, aes(x=FERx, y=FERy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("FER (scaled to year average)")+
  ylab("Effect of FER differences on alpha-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")

# # # # #
# 1.1.b.b - plants time ----
#
p1<- ggplot(splines_alpha_plant_time_sc, aes(x=LUIx, y=LUIy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("LUI (scaled to EP average)")+
  ylab("Effect of LUI differences on alpha-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p1

p2<- ggplot(splines_alpha_plant_time_sc, aes(x=MOWx, y=MOWy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("MOW (scaled to EP average)")+
  ylab("Effect of MOW differences on alpha-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_alpha_plant_time_sc, aes(x=GRAx, y=GRAy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("GRA (scaled to EP average)")+
  ylab("Effect of GRA differences on alpha-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p3


p4<- ggplot(splines_alpha_plant_time_sc, aes(x=FERx, y=FERy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("FER (scaled to EP average)")+
  ylab("Effect of FER differences on alpha-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")


# # # # #
# 1.1.b.c - herbivores space ----
#
#herbivores
#space
p1<- ggplot(splines_alpha_herb_sc_space, aes(x=LUIx, y=LUIy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("LUI (scaled to year average)")+
  ylab("Effect of LUI differences on alpha-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p1

p2<- ggplot(splines_alpha_herb_sc_space, aes(x=MOWx, y=MOWy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("MOW (scaled to year average)")+
  ylab("Effect of MOW differences on alpha-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_alpha_herb_sc_space, aes(x=GRAx, y=GRAy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("GRA (scaled to year average)")+
  ylab("Effect of GRA differences on alpha-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p3


p4<- ggplot(splines_alpha_herb_sc_space, aes(x=FERx, y=FERy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("FER (scaled to year average)")+
  ylab("Effect of FER differences on alpha-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")

# # # # #
# 1.1.b.d - herbivores time ----
# time
p1<- ggplot(splines_alpha_herb_sc_time, aes(x=LUIx, y=LUIy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("LUI (scaled to EP average)")+
  ylab("Effect of LUI differences on alpha-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p1

p2<- ggplot(splines_alpha_herb_sc_time, aes(x=MOWx, y=MOWy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("MOW (scaled to EP average)")+
  ylab("Effect of MOW differences on alpha-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_alpha_herb_sc_time, aes(x=GRAx, y=GRAy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("GRA (scaled to EP average)")+
  ylab("Effect of GRA differences on alpha-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p3


p4<- ggplot(splines_alpha_herb_sc_time, aes(x=FERx, y=FERy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("FER (scaled to EP average)")+
  ylab("Effect of FER differences on alpha-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")

# # # # #
# 1.1.b.e - sec. consumers space ----
#
# sec. consumers
# space
p1<- ggplot(splines_alpha_pred_sc_space, aes(x=LUIx, y=LUIy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("LUI (scaled to year average)")+
  ylab("Effect of LUI differences on alpha-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p1

p2<- ggplot(splines_alpha_pred_sc_space, aes(x=MOWx, y=MOWy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("MOW (scaled to year average)")+
  ylab("Effect of MOW differences on alpha-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_alpha_pred_sc_space, aes(x=GRAx, y=GRAy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("GRA (scaled to year average)")+
  ylab("Effect of GRA differences on alpha-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p3


p4<- ggplot(splines_alpha_pred_sc_space, aes(x=FERx, y=FERy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("FER (scaled to year average)")+
  ylab("Effect of FER differences on alpha-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")

# # # # #
# 1.1.b.f - sec. consumers time ----
#
#time
p1<- ggplot(splines_alpha_pred_sc_time, aes(x=LUIx, y=LUIy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("LUI (scaled to EP average)")+
  ylab("Effect of LUI differences on alpha-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p1

p2<- ggplot(splines_alpha_pred_sc_time, aes(x=MOWx, y=MOWy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("MOW (scaled to EP average)")+
  ylab("Effect of MOW differences on alpha-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_alpha_pred_sc_time, aes(x=GRAx, y=GRAy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("GRA (scaled to EP average)")+
  ylab("Effect of GRA differences on alpha-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p3


p4<- ggplot(splines_alpha_pred_sc_time, aes(x=FERx, y=FERy, colour=alpha_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 0.5), expand = c(0, 0))+
  scale_colour_hue(name="Type of Alpha Diversity", l=50)+
  xlab("FER (scaled to EP average)")+
  ylab("Effect of FER differences on alpha-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")


# # # # #
# 1.2 - BETADIVERSITY ----
#

# # # # #
# 1.2.a. - MAIN FIGURES  ----
#
# beta turnover
# subset data
library(dplyr)
subsp <- div_space[grep('bsim', div_space$beta_type),]
subti <- beta_time[grep('bsim', beta_time$beta_type),]

subti$organism<- factor(subti$organism, levels = c("plants", "herbivores", "sec.consumers"))
p2<- ggplot(subti, aes(x=LUIx, y=LUIy, colour=organism))+
  geom_line(size=1.5)+
  xlim(-2.5, 2.5)+
  scale_y_sqrt(limits = c(0, 1), expand = c(0, 0))+
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))+
  xlab("LUI (scaled to EP average)")+
  ylab("Effect of LUI differences on beta-diversity")+
  ggtitle("b - temporal land use intensity (LUI)\neffects on beta diversity")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  
  theme_classic()+
  theme(text=element_text(size = 20))

p2

subsp$organism<- factor(subsp$organism, levels = c("plants", "herbivores", "sec.consumers"))
p1<- ggplot(subsp, aes(x=LUIx, y=LUIy, colour=organism))+
  geom_line(size=1.5)+
  xlim(-2.5, 2.5)+
  scale_y_sqrt(limits = c(0, 1), expand = c(0, 0))+
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))+
  xlab("LUI (scaled to year average)")+
  ylab("Effect of LUI differences on beta-diversity")+
  ggtitle("a - spatial land use intensity (LUI)\neffects on beta diversity")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p1


p1+p2+plot_layout(guides = "collect")


# # # # #
# 1.2.a. - SUPPLEMENTARY FIGURES  ----
#
#plotting as png with 1500*1500 pixels

# # # # #
# 1.2.b.a. - plants space ----
#
#plants space

p1<- ggplot(splines_plants_space_sc, aes(x=LUIx, y=LUIy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), 
                   l=50)+
  xlab("LUI (scaled to year average)")+
  ylab("Effect of LUI differences on beta-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p1

p2<- ggplot(splines_plants_space_sc, aes(x=MOWx, y=MOWy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("MOW (scaled to year average)")+
  ylab("Effect of MOW differences on beta-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_plants_space_sc, aes(x=GRAx, y=GRAy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("GRA (scaled to year average)")+
  ylab("Effect of GRA differences on beta-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p3


p4<- ggplot(splines_plants_space_sc, aes(x=FERx, y=FERy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("FER (scaled to year average)")+
  ylab("Effect of FER differences on beta-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")


# # # # #
# 1.2.b.b. - plants time ----
#
#time
p1<- ggplot(splines_plants_yall_sc, aes(x=LUIx, y=LUIy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("LUI (scaled to EP average)")+
  ylab("Effect of LUI differences on beta-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p1

p2<- ggplot(splines_plants_yall_sc, aes(x=MOWx, y=MOWy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("MOW (scaled to EP average)")+
  ylab("Effect of MOW differences on Beta-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_plants_yall_sc, aes(x=GRAx, y=GRAy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("GRA (scaled to EP average)")+
  ylab("Effect of GRA differences on beta-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p3


p4<- ggplot(splines_plants_yall_sc, aes(x=FERx, y=FERy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("FER (scaled to EP average)")+
  ylab("Effect of FER differences on beta-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")


# # # # #
# 1.2.b.c. - herbivores space ----
#
#herbivores space
p1<- ggplot(splines_her_space_sc, aes(x=LUIx, y=LUIy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("LUI (scaled to year average)")+
  ylab("Effect of LUI differences on beta-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p1

p2<- ggplot(splines_her_space_sc, aes(x=MOWx, y=MOWy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("MOW (scaled to year average)")+
  ylab("Effect of MOW differences on beta-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_her_space_sc, aes(x=GRAx, y=GRAy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("GRA (scaled to year average)")+
  ylab("Effect of GRA differences on beta-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p3


p4<- ggplot(splines_her_space_sc, aes(x=FERx, y=FERy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("FER (scaled to year average)")+
  ylab("Effect of FER differences on beta-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")

# # # # #
# 1.2.b.d. - herbivores time ----
#
#time
p1<- ggplot(splines_her_yall_sc, aes(x=LUIx, y=LUIy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("LUI (scaled to EP average)")+
  ylab("Effect of LUI differences on beta-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p1

p2<- ggplot(splines_her_yall_sc, aes(x=MOWx, y=MOWy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("MOW (scaled to EP average)")+
  ylab("Effect of MOW differences on beta-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_her_yall_sc, aes(x=GRAx, y=GRAy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("GRA (scaled to EP average)")+
  ylab("Effect of GRA differences on beta-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p3


p4<- ggplot(splines_her_yall_sc, aes(x=FERx, y=FERy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("FER (scaled to EP average)")+
  ylab("Effect of FER differences on beta-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")

# # # # #
# 1.2.b.e. - sec. consumers space ----
#
# sec. consumers space
p1<- ggplot(splines_pred_space_sc, aes(x=LUIx, y=LUIy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("LUI (scaled to year average)")+
  ylab("Effect of LUI differences on beta-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p1

p2<- ggplot(splines_pred_space_sc, aes(x=MOWx, y=MOWy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("MOW (scaled to year average)")+
  ylab("Effect of MOW differences on beta-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_pred_space_sc, aes(x=GRAx, y=GRAy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("GRA (scaled to year average)")+
  ylab("Effect of GRA differences on beta-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p3


p4<- ggplot(splines_pred_space_sc, aes(x=FERx, y=FERy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("FER (scaled to year average)")+
  ylab("Effect of FER differences on beta-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")

# # # # #
# 1.2.b.f. - sec. consumers time----
#
#time
p1<- ggplot(splines_pred_yall_sc, aes(x=LUIx, y=LUIy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("LUI (scaled to EP average)")+
  ylab("Effect of LUI differences on beta-diversity")+
  ggtitle("a - Land use intensity (LUI)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))

p1

p2<- ggplot(splines_pred_yall_sc, aes(x=MOWx, y=MOWy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("MOW (scaled to EP average)")+
  ylab("Effect of MOW differences on beta-diversity")+
  ggtitle("b - Mowing frequency (MOW)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p2


p3<- ggplot(splines_pred_yall_sc, aes(x=GRAx, y=GRAy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("GRA (scaled to EP average)")+
  ylab("Effect of GRA differences on beta-diversity")+
  ggtitle("c - Grazing intensity (GRA)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p3


p4<- ggplot(splines_pred_yall_sc, aes(x=FERx, y=FERy, colour=beta_type))+
  geom_line(size=1.5)+
  scale_y_sqrt(limits = c(0, 3), expand = c(0, 0))+
  scale_colour_hue(name="Type of Beta Diversity",                     
                   labels = c("bsim", "b0", "b1", "b2", "b3", "b4"), l=50)+
  xlab("FER (scaled to EP average)")+
  ylab("Effect of FER differences on beta-diversity")+
  ggtitle("d - Fertilisation intensity (FER)")+
  #geom_text(aes(x = 3.9, y = 0.45, label = "EV% = 0.1"), color="black") + 
  theme_classic()+
  theme(text=element_text(size = 20))


p4

p1 + p2 + p3 + p4 + plot_layout(guides="collect")



# # # # # # # # # # # # # #
# 2 - LINEAR MODELS Barplots              ----
#

# # # # #
# 2.1. - ALPHA ----

load("./data/InputData/EV_time_LUI_alpha.RData")
load("./data/InputData/EV_space_LUI_alpha.RData")

EV_time_alpha<- EV_time
EV_space_alpha<- EV_space
# # # # #
# 2.1.a. - alpha time ----

#EV
a0_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[1],EV_time[[2]]$mLUIpp[1],EV_time[[3]]$mLUIpp[1])
a0_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[1],EV_time[[2]]$dLUIpp[1],EV_time[[3]]$dLUIpp[1])

a1_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[2],EV_time[[2]]$mLUIpp[2],EV_time[[3]]$mLUIpp[2])
a1_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[2],EV_time[[2]]$dLUIpp[2],EV_time[[3]]$dLUIpp[2])

a2_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[3],EV_time[[2]]$mLUIpp[3],EV_time[[3]]$mLUIpp[3])
a2_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[3],EV_time[[2]]$dLUIpp[3],EV_time[[3]]$dLUIpp[3])

a3_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[4],EV_time[[2]]$mLUIpp[4],EV_time[[3]]$mLUIpp[4])
a3_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[4],EV_time[[2]]$dLUIpp[4],EV_time[[3]]$dLUIpp[4])

a4_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[5],EV_time[[2]]$mLUIpp[5],EV_time[[3]]$mLUIpp[5])
a4_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[5],EV_time[[2]]$dLUIpp[5],EV_time[[3]]$dLUIpp[5])

#coefficients
a0_time_mLUI<-c(EV_time[[1]]$cofm[1],EV_time[[2]]$cofm[1],EV_time[[3]]$cofm[1])
a0_time_dLUI<-c(EV_time[[1]]$cofd[1],EV_time[[2]]$cofd[1],EV_time[[3]]$cofd[1])

a1_time_mLUI<-c(EV_time[[1]]$cofm[2],EV_time[[2]]$cofm[2],EV_time[[3]]$cofm[2])
a1_time_dLUI<-c(EV_time[[1]]$cofd[2],EV_time[[2]]$cofd[2],EV_time[[3]]$cofd[2])

a2_time_mLUI<-c(EV_time[[1]]$cofm[3],EV_time[[2]]$cofm[3],EV_time[[3]]$cofm[3])
a2_time_dLUI<-c(EV_time[[1]]$cofd[3],EV_time[[2]]$cofd[3],EV_time[[3]]$cofd[3])

a3_time_mLUI<-c(EV_time[[1]]$cofm[4],EV_time[[2]]$cofm[4],EV_time[[3]]$cofm[4])
a3_time_dLUI<-c(EV_time[[1]]$cofd[4],EV_time[[2]]$cofd[4],EV_time[[3]]$cofd[4])

a4_time_mLUI<-c(EV_time[[1]]$cofm[5],EV_time[[2]]$cofm[5],EV_time[[3]]$cofm[5])
a4_time_dLUI<-c(EV_time[[1]]$cofd[5],EV_time[[2]]$cofd[5],EV_time[[3]]$cofd[5])

# # # # #
# 2.1.b. - alpha space ----

# alpha
#
#EV
a0_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[1],EV_space[[2]]$mLUIpp[1],EV_space[[3]]$mLUIpp[1])
a0_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[1],EV_space[[2]]$dLUIpp[1],EV_space[[3]]$dLUIpp[1])

a1_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[2],EV_space[[2]]$mLUIpp[2],EV_space[[3]]$mLUIpp[2])
a1_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[2],EV_space[[2]]$dLUIpp[2],EV_space[[3]]$dLUIpp[2])

a2_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[3],EV_space[[2]]$mLUIpp[3],EV_space[[3]]$mLUIpp[3])
a2_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[3],EV_space[[2]]$dLUIpp[3],EV_space[[3]]$dLUIpp[3])

a3_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[4],EV_space[[2]]$mLUIpp[4],EV_space[[3]]$mLUIpp[4])
a3_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[4],EV_space[[2]]$dLUIpp[4],EV_space[[3]]$dLUIpp[4])

a4_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[5],EV_space[[2]]$mLUIpp[5],EV_space[[3]]$mLUIpp[5])
a4_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[5],EV_space[[2]]$dLUIpp[5],EV_space[[3]]$dLUIpp[5])

#coefficients
a0_space_mLUI<-c(EV_space[[1]]$cofm[1],EV_space[[2]]$cofm[1],EV_space[[3]]$cofm[1])
a0_space_dLUI<-c(EV_space[[1]]$cofd[1],EV_space[[2]]$cofd[1],EV_space[[3]]$cofd[1])

a1_space_mLUI<-c(EV_space[[1]]$cofm[2],EV_space[[2]]$cofm[2],EV_space[[3]]$cofm[2])
a1_space_dLUI<-c(EV_space[[1]]$cofd[2],EV_space[[2]]$cofd[2],EV_space[[3]]$cofd[2])

a2_space_mLUI<-c(EV_space[[1]]$cofm[3],EV_space[[2]]$cofm[3],EV_space[[3]]$cofm[3])
a2_space_dLUI<-c(EV_space[[1]]$cofd[3],EV_space[[2]]$cofd[3],EV_space[[3]]$cofd[3])

a3_space_mLUI<-c(EV_space[[1]]$cofm[4],EV_space[[2]]$cofm[4],EV_space[[3]]$cofm[4])
a3_space_dLUI<-c(EV_space[[1]]$cofd[4],EV_space[[2]]$cofd[4],EV_space[[3]]$cofd[4])

a4_space_mLUI<-c(EV_space[[1]]$cofm[5],EV_space[[2]]$cofm[5],EV_space[[3]]$cofm[5])
a4_space_dLUI<-c(EV_space[[1]]$cofd[5],EV_space[[2]]$cofd[5],EV_space[[3]]$cofd[5])

tlnames <- c("plants","herbivores","secondary consumers")

cols <- c("#1b9e77", "#d95f02", "#7570b3")
letmat <- matrix(letters[1:9], nrow=3,byrow=T)

#remove files as names overlaps with next files
remove(EV_time_alpha, EV_space_alpha)

# # # # #
# 2.2. - a0 MAIN FIGURE ----
#
#a0 - **main figure**


# # # # #
# 2.2.a. - a0 time ----
#time
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a0_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.07, 0.03),
             main="mean LUI")
text(x=0.02, y=y, cex=2.5,
     labels= paste(round(rev(abs(a0_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a0_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.07,0.03),
             main="delta LUI")

text(x=-0.06, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a0_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.11,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (richness, q = 0)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.b. - a0 space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a0_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.2),
             main="mean LUI")
text(x=0.1, y=y, cex=2.5,
     labels= paste(round(rev(abs(a0_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a0_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.2),
             main="delta LUI")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a0_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.6,-1.0,expression(paste("Effect on spatial differences in ",
                                alpha, " diversity (richness, q = 0)",sep="")),
     xpd=NA,cex=2.5)




# # # # #
# 2.3. - a1-4 SUPPL FIGURE ----
#
# a1-4 - **supplementary figures**
# each time and space

# # # # #
# 2.3.a. - a1 time ----
#
# time
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a1_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.07, 0.03),
             main="mean LUI")
text(x=0.02, y=y, cex=2.5,
     labels= paste(round(rev(abs(a1_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a1_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.07,0.03),
             main="delta LUI")

text(x=-0.06, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a1_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.10,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 1)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.3.b. - a1 space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a1_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.2,0.4),
             main="mean LUI")
text(x=0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a1_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a1_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.2,0.4),
             main="delta LUI")
text(x=-0.1, y=y, cex=2.5,
     labels= paste(round(rev(abs(a1_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.5,-1.0,expression(paste("Effect on spatial differences in ",
                                alpha, " diversity (q = 1)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.3.c. - a2 time ----
#
#time
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a2_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.07, 0.03),
             main="mean LUI")
text(x=0.02, y=y, cex=2.5,
     labels= paste(round(rev(abs(a2_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a2_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.07,0.03),
             main="delta LUI")

text(x=-0.06, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a2_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.10,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 2)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.3.d. - a2 space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a2_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.4,0.2),
             main="mean LUI")
text(x=0.1, y=y, cex=2.5,
     labels= paste(round(rev(abs(a2_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a2_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.2),
             main="delta LUI")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a2_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.5,-1.0,expression(paste("Effect on spatial differences in ",
                                alpha, " diversity (q = 2)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.3.e. - a3 time ----
#
# a3
#time
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a3_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T, names=rev(tlnames),
             xlim=c(-0.07, 0.03),
             main="mean LUI")
text(x=0.02, y=y, cex=2.5,
     labels= paste(round(rev(abs(a3_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a3_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.07,0.03),
             main="delta LUI")

text(x=-0.06, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a3_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.10,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 3)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.3.f.  -  a3 space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a3_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.2),
             main="mean LUI")
text(x=0.1, y=y, cex=2.5,
     labels= paste(round(rev(abs(a3_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a3_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.2),
             main="delta LUI")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a3_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.5,-1.0,expression(paste("Effect on spatial differences in ",
                                alpha, " diversity (q = 3)",sep="")),
     xpd=NA,cex=2.5)


# # # # #
# 2.3.g. - a4 time ----
#
# a4
#time
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(8,6,3,0))

y <- barplot(rev(a4_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.07, 0.03),
             main="mean LUI")
text(x=0.02, y=y, cex=2.5,
     labels= paste(round(rev(abs(a4_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(a4_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5, 
             axisnames = T,
             xlim=c(-0.07,0.03),
             main="delta LUI")

text(x=-0.06, y=z,cex=2.5, 
     labels= paste(round(rev(abs(a4_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.10,-1, expression(paste("Effect on temporal differences in ", 
                                alpha," diversity (q = 4)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.3.h. - a4 space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(8,6,3,0))


y <- barplot(rev(a4_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5, cex.axis = 2.5, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.2),
             main="mean LUI")
text(x=0.1, y=y, cex=2.5,
     labels= paste(round(rev(abs(a4_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(a4_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2.5, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.2),
             main="delta LUI")
text(x=-0.3, y=y, cex=2.5,
     labels= paste(round(rev(abs(a4_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.5,-1.0,expression(paste("Effect on spatial differences in ",
                                alpha, " diversity (q = 4)",sep="")),
     xpd=NA,cex=2.5)




# # # # #
# 2.2. - BETA ----
#
load("./data/InputData/EV_time_LUI_beta.RData")
load("./data/InputData/EV_space_LUI_beta.RData")

#time
#EV
bsim_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[1],EV_time[[2]]$mLUIpp[1],EV_time[[3]]$mLUIpp[1])
bsim_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[1],EV_time[[2]]$dLUIpp[1],EV_time[[3]]$dLUIpp[1])

b0_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[2],EV_time[[2]]$mLUIpp[2],EV_time[[3]]$mLUIpp[2])
b0_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[2],EV_time[[2]]$dLUIpp[2],EV_time[[3]]$dLUIpp[2])

b1_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[3],EV_time[[2]]$mLUIpp[3],EV_time[[3]]$mLUIpp[3])
b1_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[3],EV_time[[2]]$dLUIpp[3],EV_time[[3]]$dLUIpp[3])

b2_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[4],EV_time[[2]]$mLUIpp[4],EV_time[[3]]$mLUIpp[4])
b2_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[4],EV_time[[2]]$dLUIpp[4],EV_time[[3]]$dLUIpp[4])

b3_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[5],EV_time[[2]]$mLUIpp[5],EV_time[[3]]$mLUIpp[5])
b3_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[5],EV_time[[2]]$dLUIpp[5],EV_time[[3]]$dLUIpp[5])

b4_time_mLUI_EV<-c(EV_time[[1]]$mLUIpp[6],EV_time[[2]]$mLUIpp[6],EV_time[[3]]$mLUIpp[6])
b4_time_dLUI_EV<-c(EV_time[[1]]$dLUIpp[6],EV_time[[2]]$dLUIpp[6],EV_time[[3]]$dLUIpp[6])

#coefficients
bsim_time_mLUI<-c(EV_time[[1]]$cofm[1],EV_time[[2]]$cofm[1],EV_time[[3]]$cofm[1])
bsim_time_dLUI<-c(EV_time[[1]]$cofd[1],EV_time[[2]]$cofd[1],EV_time[[3]]$cofd[1])

b0_time_mLUI<-c(EV_time[[1]]$cofm[2],EV_time[[2]]$cofm[2],EV_time[[3]]$cofm[2])
b0_time_dLUI<-c(EV_time[[1]]$cofd[2],EV_time[[2]]$cofd[2],EV_time[[3]]$cofd[2])

b1_time_mLUI<-c(EV_time[[1]]$cofm[3],EV_time[[2]]$cofm[3],EV_time[[3]]$cofm[3])
b1_time_dLUI<-c(EV_time[[1]]$cofd[3],EV_time[[2]]$cofd[3],EV_time[[3]]$cofd[3])

b2_time_mLUI<-c(EV_time[[1]]$cofm[4],EV_time[[2]]$cofm[4],EV_time[[3]]$cofm[4])
b2_time_dLUI<-c(EV_time[[1]]$cofd[4],EV_time[[2]]$cofd[4],EV_time[[3]]$cofd[4])

b3_time_mLUI<-c(EV_time[[1]]$cofm[5],EV_time[[2]]$cofm[5],EV_time[[3]]$cofm[5])
b3_time_dLUI<-c(EV_time[[1]]$cofd[5],EV_time[[2]]$cofd[5],EV_time[[3]]$cofd[5])

b4_time_mLUI<-c(EV_time[[1]]$cofm[6],EV_time[[2]]$cofm[6],EV_time[[3]]$cofm[6])
b4_time_dLUI<-c(EV_time[[1]]$cofd[6],EV_time[[2]]$cofd[6],EV_time[[3]]$cofd[6])


#space
#EV
bsim_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[1],EV_space[[2]]$mLUIpp[1],EV_space[[3]]$mLUIpp[1])
bsim_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[1],EV_space[[2]]$dLUIpp[1],EV_space[[3]]$dLUIpp[1])

b0_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[2],EV_space[[2]]$mLUIpp[2],EV_space[[3]]$mLUIpp[2])
b0_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[2],EV_space[[2]]$dLUIpp[2],EV_space[[3]]$dLUIpp[2])

b1_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[3],EV_space[[2]]$mLUIpp[3],EV_space[[3]]$mLUIpp[3])
b1_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[3],EV_space[[2]]$dLUIpp[3],EV_space[[3]]$dLUIpp[3])

b2_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[4],EV_space[[2]]$mLUIpp[4],EV_space[[3]]$mLUIpp[4])
b2_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[4],EV_space[[2]]$dLUIpp[4],EV_space[[3]]$dLUIpp[4])

b3_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[5],EV_space[[2]]$mLUIpp[5],EV_space[[3]]$mLUIpp[5])
b3_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[5],EV_space[[2]]$dLUIpp[5],EV_space[[3]]$dLUIpp[5])

b4_space_mLUI_EV<-c(EV_space[[1]]$mLUIpp[6],EV_space[[2]]$mLUIpp[6],EV_space[[3]]$mLUIpp[6])
b4_space_dLUI_EV<-c(EV_space[[1]]$dLUIpp[6],EV_space[[2]]$dLUIpp[6],EV_space[[3]]$dLUIpp[6])

#coefficients
bsim_space_mLUI<-c(EV_space[[1]]$cofm[1],EV_space[[2]]$cofm[1],EV_space[[3]]$cofm[1])
bsim_space_dLUI<-c(EV_space[[1]]$cofd[1],EV_space[[2]]$cofd[1],EV_space[[3]]$cofd[1])

b0_space_mLUI<-c(EV_space[[1]]$cofm[2],EV_space[[2]]$cofm[2],EV_space[[3]]$cofm[2])
b0_space_dLUI<-c(EV_space[[1]]$cofd[2],EV_space[[2]]$cofd[2],EV_space[[3]]$cofd[2])

b1_space_mLUI<-c(EV_space[[1]]$cofm[3],EV_space[[2]]$cofm[3],EV_space[[3]]$cofm[3])
b1_space_dLUI<-c(EV_space[[1]]$cofd[3],EV_space[[2]]$cofd[3],EV_space[[3]]$cofd[3])

b2_space_mLUI<-c(EV_space[[1]]$cofm[4],EV_space[[2]]$cofm[4],EV_space[[3]]$cofm[4])
b2_space_dLUI<-c(EV_space[[1]]$cofd[4],EV_space[[2]]$cofd[4],EV_space[[3]]$cofd[4])

b3_space_mLUI<-c(EV_space[[1]]$cofm[5],EV_space[[2]]$cofm[5],EV_space[[3]]$cofm[5])
b3_space_dLUI<-c(EV_space[[1]]$cofd[5],EV_space[[2]]$cofd[5],EV_space[[3]]$cofd[5])

b4_space_mLUI<-c(EV_space[[1]]$cofm[6],EV_space[[2]]$cofm[6],EV_space[[3]]$cofm[6])
b4_space_dLUI<-c(EV_space[[1]]$cofd[6],EV_space[[2]]$cofd[6],EV_space[[3]]$cofd[6])


tlnames <- c("plants","herbivores","secondary consumers")

cols <- c("#1b9e77", "#d95f02", "#7570b3")
letmat <- matrix(letters[1:9], nrow=3,byrow=T)

#remove files as they overlap with next ones in their name
remove(EV_time, EV_space)

# # # # #
# 2.2.a. - turnover time ----
#
#Turnover
#time
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(bsim_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.04, 0.06),
             main="mean LUI")
text(x=0.05, y=y,cex=2.5,
     labels= paste(round(rev(abs(bsim_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(bsim_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.04,0.06),
             main="delta LUI")

text(x=-0.03, y=z, cex=2.5,
     labels= paste(round(rev(abs(bsim_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.07,-1,expression(paste("Effect on temporal ",beta, " diversity (turnover)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.b. - turnover space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(bsim_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.3),
             main="mean LUI")
text(x=0.2, y=y, cex=2.5,
     labels= paste(round(rev(abs(bsim_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(bsim_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.3),
             main="delta LUI")
text(x=-0.5, y=y, cex=2.5,
     labels= paste(round(rev(abs(bsim_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.6,-1.0,
     expression(paste("Effect on spatial ",beta, " diversity (turnover)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.c. - b0 time ----
#
#b0
#time
dev.off()
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b0_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.04, 0.06),
             main="mean LUI")
text(x=0.05, y=y,cex=2.5,
     labels= paste(round(rev(abs(b0_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(b0_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.04,0.06),
             main="delta LUI")

text(x=-0.03, y=z, cex=2.5,
     labels= paste(round(rev(abs(b0_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.07,-1,expression(paste("Effect on temporal ",beta, " diversity (q=0)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.d. - b0 space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b0_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.3),
             main="mean LUI")
text(x=0.2, y=y, cex=2.5,
     labels= paste(round(rev(abs(b0_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b0_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.3),
             main="delta LUI")
text(x=-0.5, y=y, cex=2.5,
     labels= paste(round(rev(abs(b0_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.6,-1.0,
     expression(paste("Effect on spatial ",beta, " diversity (q=0)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.e. - b1 time ----
#
#b1
#time
dev.off()
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b1_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.04, 0.06),
             main="mean LUI")
text(x=0.05, y=y,cex=2.5,
     labels= paste(round(rev(abs(b1_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(b1_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.04,0.06),
             main="delta LUI")

text(x=-0.03, y=z, cex=2.5,
     labels= paste(round(rev(abs(b1_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.07,-1,expression(paste("Effect on temporal ",beta, " diversity (q=1)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.f. - b1 space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b1_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.3),
             main="mean LUI")
text(x=0.2, y=y, cex=2.5,
     labels= paste(round(rev(abs(b1_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b1_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.3),
             main="delta LUI")
text(x=-0.5, y=y, cex=2.5,
     labels= paste(round(rev(abs(b1_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.6,-1.0,
     expression(paste("Effect on spatial ",beta, " diversity (q=1)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.g. - b2 time ----
#
#b2
#time
dev.off()
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b2_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.04, 0.06),
             main="mean LUI")
text(x=0.05, y=y,cex=2.5,
     labels= paste(round(rev(abs(b2_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(b2_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.04,0.06),
             main="delta LUI")

text(x=-0.03, y=z, cex=2.5,
     labels= paste(round(rev(abs(b2_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.07,-1,expression(paste("Effect on temporal ",beta, " diversity (q=2)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.h. - b2 space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b2_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.3),
             main="mean LUI")
text(x=0.2, y=y, cex=2.5,
     labels= paste(round(rev(abs(b2_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b2_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.3),
             main="delta LUI")
text(x=-0.5, y=y, cex=2.5,
     labels= paste(round(rev(abs(b2_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.6,-1.0,
     expression(paste("Effect on spatial ",beta, " diversity (q=2)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.i. - b3 time ----
#
#b3
#time
dev.off()
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b3_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.04, 0.06),
             main="mean LUI")
text(x=0.05, y=y,cex=2.5,
     labels= paste(round(rev(abs(b3_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(b3_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.04,0.06),
             main="delta LUI")

text(x=-0.03, y=z, cex=2.5,
     labels= paste(round(rev(abs(b3_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.07,-1,expression(paste("Effect on temporal ",beta, " diversity (q=3)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.j. - b3 space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b3_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.3),
             main="mean LUI")
text(x=0.2, y=y, cex=2.5,
     labels= paste(round(rev(abs(b3_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b3_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.3),
             main="delta LUI")
text(x=-0.5, y=y, cex=2.5,
     labels= paste(round(rev(abs(b3_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.6,-1.0,
     expression(paste("Effect on spatial ",beta, " diversity (q=3)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.k. - b4 time ----
#
#b4
#time
dev.off()
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b4_time_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.04, 0.06),
             main="mean LUI")
text(x=0.05, y=y,cex=2.5,
     labels= paste(round(rev(abs(b4_time_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

z <- barplot(rev(b4_time_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.04,0.06),
             main="delta LUI")

text(x=-0.03, y=z, cex=2.5,
     labels= paste(round(rev(abs(b4_time_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.07,-1,expression(paste("Effect on temporal ",beta, " diversity (q=4)",sep="")),
     xpd=NA,cex=2.5)

# # # # #
# 2.2.l. - b4 space ----
#
#space
par(mfrow=c(1,2))
par(mar=c(2,15,1,2))
par(oma=c(6,8,3,0))

y <- barplot(rev(b4_space_mLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,names=rev(tlnames),
             xlim=c(-0.4,0.3),
             main="mean LUI")
text(x=0.2, y=y, cex=2.5,
     labels= paste(round(rev(abs(b4_space_mLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))


z <- barplot(rev(b4_space_dLUI), col = rev(cols),las=1,
             horiz=T, cex.main=2.5,cex.axis = 2, cex.names=2.5,
             axisnames = T,
             xlim=c(-0.4,0.3),
             main="delta LUI")
text(x=-0.5, y=y, cex=2.5,
     labels= paste(round(rev(abs(b4_space_dLUI_EV)), 0), "%", sep=""),
     adj=c(0.5,0.5))

## combined x axis
text(-0.6,-1.0,
     expression(paste("Effect on spatial ",beta, " diversity (q=4)",sep="")),
     xpd=NA,cex=2.5)


# these plots above can be replicated for the other LUI components

# # # # # # # # # # # # # #
# 3 - SCATTERPLOTS              ----
#
# Overview scatterplot - average dissimilarity in space and time
# idea: #scatterplot group means - space vs time
#
# # # # #
# 3.x.x. - create all_region_means.RData ----
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # #code to produce all_region_means.RData # # # # # # # # # # #
# load packages
library(data.table)

#TODO HERE

# set working directory to 4. plotting results
# load required datasets as produced/loaded above (only if you start running the script here)

load("div_temp_fig.RData")
load("div_space_fig.RData")

load("./data/InputData/pwise_time_plants.RData")
load("./data/InputData/pwise_time_herb.RData")
load("./data/InputData/pwise_time_pred.RData")

load("./data/InputData/pwise_space_plants.RData")
load("./data/InputData/pwise_space_herb.RData")
load("./data/InputData/pwise_space_pred.RData")

# set working directory to R Functions
# and the we link the functions to plott results to this script
source("./functions_plotting_results.R")

# # # # dataset with all beta measures, averaged per region
## spatial data only within region averages

#TODO remove function, keep explanation
# calc_se <- function(x, N) sqrt(var(x, na.rm = T) / N) # function to calc standard error
# need to input N : the number of unique observations (number of plots * number of year)
div_temp$N_unique_obs <- length(unique(pwise_time_insects$EPy1))
div_temp$N_unique_obs[which(div_temp$group == "plant")] <- length(unique(pwise_time_plants$EPy1))
div_space$N_unique_obs <- length(unique(pwise_space_insects$EPy1))
div_space$N_unique_obs[which(div_space$group == "plant")] <- length(unique(pwise_space_plants$EPy1))

# create spatial dataset with only-within-region
div_space_within_reg <- div_space[which(div_space$reg %in% c("AEG.AEG", "HEG.HEG", "SEG.SEG")), ]
div_space_within_reg$reg <- sapply(strsplit(as.character(div_space_within_reg$reg), "\\."), "[", 2) # rename reg to match temporal dataset
rm(div_space)

# functions use data.table solution --> therefore need to change class.
div_temp <- data.table(div_temp)
div_space_within_reg <- data.table(div_space_within_reg)


# CREATE MEAN AND STANDARD ERROR PER REGION PER GROUP ----
# BSIM
temp <- calc_temporal_regionwise_mean_sd("Bsim_time", div_temp, "bsim")
bsim <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "Bsim_space", div_space_within_reg, "bsim")

#B0
temp <- calc_temporal_regionwise_mean_sd("B0_time", div_temp, "b0")
b0 <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "B0_space", div_space_within_reg, "b0")

#B1
temp <- calc_temporal_regionwise_mean_sd("B1_time", div_temp, "b1")
b1 <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "B1_space", div_space_within_reg, "b1")

#B1
temp <- calc_temporal_regionwise_mean_sd("B2_time", div_temp, "b2")
b2 <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "B2_space", div_space_within_reg, "b2")

#b3
temp <- calc_temporal_regionwise_mean_sd("B3_time", div_temp, "b3")
b3 <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "B3_space", div_space_within_reg, "b3")

#b4
temp <- calc_temporal_regionwise_mean_sd("B4_time", div_temp, "b4")
b4 <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "B4_space", div_space_within_reg, "b4")

#a0
temp <- calc_temporal_regionwise_mean_sd("A0_time", div_temp, "a0")
a0 <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "A0_space", div_space_within_reg, "a0")

#a1
temp <- calc_temporal_regionwise_mean_sd("A1_time", div_temp, "a1")
a1 <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "A1_space", div_space_within_reg, "a1")

#a2
temp <- calc_temporal_regionwise_mean_sd("A2_time", div_temp, "a2")
a2 <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "A2_space", div_space_within_reg, "a2")

#a3
temp <- calc_temporal_regionwise_mean_sd("A3_time", div_temp, "a3")
a3 <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "A3_space", div_space_within_reg, "a3")

#a4
temp <- calc_temporal_regionwise_mean_sd("A4_time", div_temp, "a4")
a4 <- calc_spatial_regionwise_mean_sd_and_add_to_temporal(temp, "A4_space", div_space_within_reg, "a4")


reg<- as.factor(rep(c("AEG", "HEG", "SEG"), each=3))

# # # # # # # # # # # # # #
# 3 - SCATTERPLOTS              ----
#
library(ggplot2)
library(cowplot)

# set working directory to R Functions
# and the we link the functions to plott results to this script
source("./functions_plotting_results.R")

# always : plants in green, herbivores in orange and secondary consumers in violet
#  note : less error-prone to include colors of plotting in table
coltable <- data.table(group = c("plant", "herbivore", "secondary consumer"), col_by_group = c("#1b9e77", "#d95f02", "#7570b3"))
for(d in c(paste("b", c("sim", seq(0, 4)), sep = ""), paste("a", seq(0, 4), sep = ""))){
  # for loop performs below code for all datasets
  # bsim <- merge(bsim, coltable, by = c("group")) # add colors to bsim
  temp <- merge(get(d), coltable, by = c("group"))
  assign(d, temp)
}
# plotting parameters
errorbar_width <- 0.01

# plotting
p_b0 <- generate_time_space_comparison_scatterplot(dat = b0, dat_name = "b0", 
                                                   pretty_name = "q=0")
p_b1 <- generate_time_space_comparison_scatterplot(dat = b1, dat_name = "b1", 
                                                   pretty_name = "q=1")
p_b2 <- generate_time_space_comparison_scatterplot(dat = b2, dat_name = "b2",
                                                   pretty_name = "q=2")
p_b3 <- generate_time_space_comparison_scatterplot(dat = b2, dat_name = "b3",
                                                   pretty_name = "q=3")
p_b4 <- generate_time_space_comparison_scatterplot(dat = b2, dat_name = "b4",
                                                   pretty_name = "q=4")

p_a0 <- generate_time_space_comparison_scatterplot(dat = b0, dat_name = "a0", 
                                                   pretty_name = "species richness, q=0")
p_a1 <- generate_time_space_comparison_scatterplot(dat = b1, dat_name = "a1", 
                                                   pretty_name = "q=1")
p_a2 <- generate_time_space_comparison_scatterplot(dat = b2, dat_name = "a2",
                                                   pretty_name = "q=2")
p_a3 <- generate_time_space_comparison_scatterplot(dat = b2, dat_name = "a3",
                                                   pretty_name = "q=3")
p_a4 <- generate_time_space_comparison_scatterplot(dat = b2, dat_name = "a4",
                                                   pretty_name = "q=4")
                                                   
# if you use cowplot, you can elegantly combine several plots
cowplot::plot_grid(p_a0, p_bsim,  labels = c("(a)", "(b)"))

plot_grid(p_bsim, p_b0, p_b1, p_b2, labels = c("(a)", "(b)", "(c)", "(d)"))
# see this tutorial (I use it a lot): https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html


