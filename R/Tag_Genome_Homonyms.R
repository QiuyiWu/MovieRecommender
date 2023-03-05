#some analyis here
tag = read.csv("../../Data/genome-tags.csv")
tags = tag$tag
#1084 is "violence", 1085 is "violent"
homonyms = NULL
homonyms[[1]] = c(1,2) #007 and 007 (series) ? 
homonyms[[2]] = c(9,13) # 1980's and 80's
homonyms[[3]] = c(14,15) #aardman and aardman animations
homonyms[[4]] = c(19,20) #action and action-packed
homonyms[[5]] = c(52,53) #alteranate rality and alternate univesre
homonyms[[6]] = c(59,60) #"android(s)/cyborg(s)" "androids"    
homonyms[[7]] = c(81,82,83) #assasin, assassination, assasisns
homonyms[[8]] = c(87,88) #australia and australian
homonyms[[9]] = c(133,134,135) #biographical, biography, biopic
homonyms[[10]] = c(145,146) #blood and bloody
homonyms[[11]] = c(155,156) #boring and boring!
homonyms[[12]] = c(168,169) #brutal and brutality
homonyms[[13]] = c(208,209) #christian and christianity
homonyms[[14]] = c(221,222) #clones and cloning
homonyms[[15]] = c(224,225) #coen bros and coen brothers
homonyms[[16]] = c(231,232,233,234) #comic, comic book, comic book adaption, comics
homonyms[[17]] = c(235,236) #coming of age and coming-of-age
homonyms[[18]] = c(247,248) #con artists and con men
homonyms[[19]] = c(261,262,263) #court, courtroom, courtroom drama
homonyms[[20]] = c(283,284) #dance, dancing
homonyms[[21]] = c(321,322) #dragon, dragons
homonyms[[22]] = c(329,330) #drug abuse, drug addictoin
homonyms[[23]] = c(357,358) #environment, environmental
homonyms[[24]] = c(371,372) #fairy tale, fairy tales
homonyms[[25]] = c(384,385) #father son relationship, father-son-relationship
homonyms[[26]] = c(387,388) #fell good movie, feel-good
homonyms[[27]] = c(414,415) #fun, fun movie
homonyms[[28]] = c(422,423,424) #gangs, gangster, gangsters
homonyms[[29]] = c(425,426) #gay, gay character
homonyms[[30]] = c(427,428) #geek, geeks
homonyms[[31]] = c(435,436) #ghost, ghosts/afterlife
homonyms[[32]] = c(456,457,458) #gore, goretatsic, gory
homonyms[[33]] = c(484,485) #hackers, hacking
homonyms[[34]] = c(505,506) #hilarious, hillarious 
homonyms[[35]] = c(508,509) #historical, history
homonyms[[36]] = c(528,529) #humor, humorous
homonyms[[37]] = c(548,549) #inspirational, inspiring
homonyms[[38]] = c(597,598) #lawyer, lawyers
homonyms[[39]] = c(635,636) #math, mathematics
homonyms[[40]] = c(663,664) #monster, monsters
homonyms[[41]] = c(696,697) #nazi, nazis
homonyms[[42]] = c(714,715) #nonlinear, non-linear
homonyms[[43]] = c(766,767) #paranoia, paranoid
homonyms[[44]] = c(785,786) #pixar, pixar animation
homonyms[[45]] = c(802,803) #postapocalyptic and post apocalypitc
homonyms[[46]] = c(848,849) #remade and remake
homonyms[[47]] = c(860,861) #robot and robots
homonyms[[48]] = c(863,864) #romance, romantic
homonyms[[49]] = c(877,878) #satire, satirical
homonyms[[50]] = c(886,887,889,890) #sci fi, sci-fi, science fiction, scifi
homonyms[[51]] = c(929,930) #slow, slow-paced
homonyms[[52]] = c(955,960,961) #spies, spy, spying
homonyms[[53]] = c(969,970) #stop motion, stop-motion
homonyms[[54]] = c(987,988,989,990) #super hero, super-hero, superhero, superheroes
homonyms[[55]] = c(999,1000) #suspense, suspenseful
homonyms[[56]] = c(1005,1006) #sword fight, sword-fighting
homonyms[[57]] = c(1013, 1015, 1016,1017) #teen, teenagers, teenager, teens
homonyms[[58]] = c(1067,1069) #vampire, vampires
homonyms[[59]] = c(1074,1075,1077,1078) #video game, video game adaptation, video games, videogames
homonyms[[60]] = c(1082,1083) #vigilante, vigilantism
homonyms[[61]] = c(1084,1085) #violence, violent
homonyms[[62]] = c(1105,1106) #werefolf, werewolves
homonyms[[63]] = c(1112,1113) #witch, witches
homonyms[[64]] = c(1121,1126) #world war ii, wwii
homonyms[[65]] = c(1127,1128) #zombie, zombies

load("../data/tag_mat.rda")

cor_vec = NULL
for(i in 1:length(homonyms)){
  temp = as.vector(cor(tag_mat[,homonyms[[i]]]))
  cor_vec = c(cor_vec,temp)
  if(sum(which(temp<0.1)) > 0){
    print(i)
  }
}

cor_vec= cor_vec[which(cor_vec < 0.99999999)]
hist(cor_vec)
